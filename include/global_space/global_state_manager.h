#pragma once
#include <bitset>

#include "local_space/neutral_atom_manifold.h"

#include "global_space/rwa_frame.h"
#include "global_space/subspace_operations.h"
#include "global_space/interaction_picture.h"

#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "hamiltonian/two_atom_rydberg.h"
#include "solvers/crank_nicolson_solver.h"


namespace KetCat
{
	template <natural_t QubitCount, NeutralAtomTypeConfig Config>
	class GlobalStateManager
	{
		using ConfigType = std::remove_cvref_t<decltype(Config)>;

		/// @brief Helper for managing the global state vector and mapping between local operations and global state.
		/// As it's a static class, we only need to define the type alias here.
		using SubspaceManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;

		/// @brief Physical manifold defining the atomic structure and seed states.
		NeutralAtomManifold<Config> m_Manifold;
        
        /// @brief Type alias for the Quantum Picture converter helper class
        using PictureConverter = InteractionPictureStateTransformer<typename SubspaceManager::FullHilbertSpace>;

        /// @brief The currently active RWA frame.
        RwaFrame<ConfigType::LevelCount> m_CurrentRwaFrame;

		/// @brief Full system state vector in the product Hilbert space, in Interaction picture
        /// @details Used during the operations, updated most of the time
		StateVector<typename SubspaceManager::FullHilbertSpace, QuantumPicture::Dirac>
            m_GlobalStateVectorDirac{};

        /// @brief Full system state vector in the product Hilbert space, in Schrödinger picture
        /// @details Used during the operations, updated most of the time
        StateVector<typename SubspaceManager::FullHilbertSpace, QuantumPicture::Schrodinger>
            m_GlobalStateVectorSchrodinger{};

    public:
		GlobalStateManager(std::bitset<QubitCount> initialState)
            : m_CurrentRwaFrame(eigenenergies_t<ConfigType::LevelCount>{}, TwoPhotonDrive{})
		{
            m_GlobalStateVectorDirac =
                SubspaceManager::basisStateFromBitstring(initialState,
					ConfigType::Logical0Level, ConfigType::Logical1Level);
		}


        /// @brief Perform a single step of the Crank-Nicolson evolution.
        ///
        /// @details
        ///    1. Constructs the local RWA Hamiltonian Ĥ(t).
        ///    2. Builds the unitary propagator U(Δt) via Crank-Nicolson.
        ///    3. Applies U(Δt) to the target qubit in the global state vector.
        void evolveOneQubitGlobalState(const TwoPhotonDrive& lasers, const natural_t targetAtom)
        {
            static const eigenenergies_t<ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            static const RwaFrame<ConfigType::LevelCount>
                RwaFrame(HartreeEnergies, lasers);

            static const eigenenergies_t<ConfigType::LevelCount>
                RwaFrameEnergies = RwaFrame.generateGlobalRwaEnergies<natural_t{1}>();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(HartreeEnergies, DipoleMatrix, lasers);

            static CrankNicolsonSolver<ConfigType::LevelCount, LinearSolverBackend::ThomasTridiagonal> Solver;

            if (TimeMaster::Clock().isInstructionStart())
            {
                prepareNewInstructionState(RwaFrame);
            }

            Hamiltonian.updateMainDiagonal(lasers);
            Hamiltonian.updateOffDiagonal(lasers);
            
            auto H = InteractionPictureHamiltonian<ConfigType::LevelCount, 1U>
                ::transform(Hamiltonian.getMatrix(), RwaFrameEnergies,
                    TimeMaster::Clock().getGlobalTime());

            Solver.updateMatrices(H, TimeMaster::Clock().getTimeStep());

            std::array<natural_t, 1> targets = { targetAtom };
            SubspaceManager::template performTimeEvolution<1>(Solver, m_GlobalStateVectorDirac, targets);
            m_GlobalStateVectorDirac.m_TimeStamp = TimeMaster::Clock().getGlobalTime();
        }

        void evolveTwoQubitGlobalState(const TwoPhotonDrive& lasers, const natural_t controlAtom, const natural_t targetAtom)
            requires (QubitCount >= 2)
        {
            static const eigenenergies_t<ConfigType::LevelCount>
                HartreeEnergies = m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount>
                DipoleMatrix = m_Manifold.getDipoleMatrix();

            static const RwaFrame<ConfigType::LevelCount>
                RwaFrame(HartreeEnergies, lasers);

            static const eigenenergies_t<ConstexprMath::pow(ConfigType::LevelCount, natural_t{2})>
                RwaFrameEnergies = RwaFrame.generateGlobalRwaEnergies<natural_t{2}>();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                SingleAtomExcitation(HartreeEnergies, DipoleMatrix, lasers);

            static TwoAtomRydbergBlockade<ConfigType::LevelCount>
                RydbergBlockade(Units::MeterToAtomicLength * 1E-9,
                    ConfigType::RydbergLevel, HartreeEnergies, DipoleMatrix);

            static CrankNicolsonSolver<ConfigType::LevelCount,
                LinearSolverBackend::FiveBandGaussianElimination> Solver;

            if (TimeMaster::Clock().isInstructionStart())
            {
                prepareNewInstructionState(RwaFrame);
            }

            SingleAtomExcitation.updateMainDiagonal(lasers);
            SingleAtomExcitation.updateOffDiagonal(lasers);
            auto SingleAtomHamiltonian = SingleAtomExcitation.getMatrix(); 

            RydbergBlockade.updateMatrix(SingleAtomHamiltonian, SingleAtomHamiltonian);
            auto H = InteractionPictureHamiltonian<ConfigType::LevelCount, 2U>
                ::transform(RydbergBlockade.getMatrix(), RwaFrameEnergies,
                    TimeMaster::Clock().getGlobalTime());

            Solver.updateMatrices(H, TimeMaster::Clock().getTimeStep());

            std::array<natural_t, 2> targets = { controlAtom, targetAtom };
            SubspaceManager::template performTimeEvolution<2>(Solver, m_GlobalStateVectorDirac, targets);
            m_GlobalStateVectorDirac.m_TimeStamp = TimeMaster::Clock().getGlobalTime();
        }

        const auto& getStateVector()
        {
            if (m_GlobalStateVectorDirac.m_TimeStamp > m_GlobalStateVectorSchrodinger.m_TimeStamp)
            {
                const eigenenergies_t Energies = m_CurrentRwaFrame.generateGlobalRwaEnergies<QubitCount>();
                m_GlobalStateVectorSchrodinger =
                    PictureConverter::toSchrodingerPicture(m_GlobalStateVectorDirac, Energies);
            }
            return m_GlobalStateVectorSchrodinger;
        }

        const auto& getManifold() const
        {
            return m_Manifold;
        }

    private:
        void prepareNewInstructionState(const decltype(m_CurrentRwaFrame)& newRwaFrame)
        {
            const eigenenergies_t CurrentEnergies = m_CurrentRwaFrame.generateGlobalRwaEnergies<QubitCount>();
            const eigenenergies_t NewEnergies = newRwaFrame.generateGlobalRwaEnergies<QubitCount>();

            m_GlobalStateVectorSchrodinger =
                PictureConverter::toSchrodingerPicture(m_GlobalStateVectorDirac, CurrentEnergies);

            m_GlobalStateVectorDirac =
                PictureConverter::toDiracPicture(m_GlobalStateVectorSchrodinger, NewEnergies);

            m_CurrentRwaFrame = newRwaFrame;
        }
	};
}