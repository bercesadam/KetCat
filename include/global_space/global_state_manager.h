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
    /// @brief Monolithic state manager tracking the joint multi-atom quantum state across frame pictures.
    /// @details Orchestrates time-dependent unitary evolution of the complete physical register. 
    /// It bridges single-qubit and multi-qubit local Hamiltonians to the global product space 
    /// H_tot = H_1 ⊗ H_2 ⊗ ... ⊗ H_n via subspace mapping utilities, 
    /// managing explicit transformations between the Dirac (interaction) and Schrödinger pictures.
    /// @tparam QubitCount Total number of physical atoms layouted in the quantum registry.
    /// @tparam Config Atomic species configuration specifying energy manifolds and electronic levels.
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

        /// @brief The RWA frame tracking atomic reference phase configurations.
        RwaFrame<ConfigType::LevelCount> m_RwaFrame;

        /// @brief Full system state vector in the product Hilbert space, in Interaction picture
        /// @details Used during the operations, updated most of the time
        StateVector<typename SubspaceManager::FullHilbertSpace, QuantumPicture::Dirac>
            m_GlobalStateVectorDirac{};

        /// @brief Full system state vector in the product Hilbert space, in Schrödinger picture
        /// @details Used during the operations, updated most of the time
        StateVector<typename SubspaceManager::FullHilbertSpace, QuantumPicture::Schrodinger>
            m_GlobalStateVectorSchrodinger{};

    public:
        /// @brief Constructs the global state manager and maps the initial binary register state.
        /// @param initialState Bitstring mapping out the active qubits to populate as |b_1,b_2,...,b_n⟩.
        GlobalStateManager(std::bitset<QubitCount> initialState)
            : m_RwaFrame(m_Manifold.getHartreeEnergies())
        {
            m_GlobalStateVectorDirac =
                SubspaceManager::basisStateFromBitstring(initialState,
                    ConfigType::Logical0Level, ConfigType::Logical1Level);

            m_GlobalStateVectorSchrodinger.m_StateVector = m_GlobalStateVectorDirac.m_StateVector;
        }


        /// @brief Perform a single step of the Crank-Nicolson evolution for a single localized target atom.
        ///
        /// @details
        ///    1. Constructs the local single-atom RWA Hamiltonian H_RWA(t).
        ///    2. Transforms H_RWA(t) to the Dirac interaction picture framework.
        ///    3. Builds the unitary propagator U(Δt) via Crank-Nicolson.
        ///    4. Embeds and applies U(Δt) on the designated target atom within the global state vector.
        ///
        /// @param lasers Physical two-photon drive field pulse configuration containing active Rabi frequencies and detunings.
        /// @param targetAtom Index of the targeted physical atom to undergo coherent excitation.
        void evolveOneQubitGlobalState(const TwoPhotonDrive& lasers, const natural_t targetAtom)
        {
            static const eigenenergies_t<ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            static const eigenenergies_t<ConfigType::LevelCount>
                RwaFrameEnergies = m_RwaFrame.getSingleRwaEnergies();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(HartreeEnergies, m_RwaFrame, DipoleMatrix);

            static CrankNicolsonSolver<ConfigType::LevelCount, LinearSolverBackend::ThomasTridiagonal> Solver;

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

        /// @brief Performs a single step of the Crank-Nicolson evolution for two coupled interacting atoms.
        ///
        /// @details
        ///    1. Constructs the single-atom excitation profiles.
        ///    2. Combines them into a 2D Hamiltonian matrix adding the non-local van der Waals penalty V_vdW.
        ///    3. Performs interaction-picture transformations using the multi-qubit joint RWA basis.
        ///    4. Numerically integrates the state vector utilizing the specialized pentadiagonal matrix solver.
        ///
        /// @param lasers Physical two-photon drive field pulse configuration containing active Rabi frequencies and detunings.
        /// @param controlAtom Index of the first interacting atom (control terminal).
        /// @param targetAtom Index of the second interacting atom (target terminal).
        void evolveTwoQubitGlobalState(const TwoPhotonDrive& lasers, const natural_t controlAtom, const natural_t targetAtom)
            requires (QubitCount >= 2)
        {
            static const eigenenergies_t<ConfigType::LevelCount>
                HartreeEnergies = m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount>
                DipoleMatrix = m_Manifold.getDipoleMatrix();

            static const eigenenergies_t < ConstexprMath::pow(ConfigType::LevelCount, natural_t{ 2 }) >
                RwaFrameEnergies = m_RwaFrame.generateGlobalRwaEnergies < natural_t{ 2 } > ();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                SingleAtomExcitation(HartreeEnergies, m_RwaFrame, DipoleMatrix);

            static TwoAtomRydbergBlockade<ConfigType::LevelCount>
                RydbergBlockade(Units::MeterToAtomicLength * 1E-8,
                    ConfigType::RydbergLevel, HartreeEnergies, DipoleMatrix);

            static CrankNicolsonSolver<ConfigType::LevelCount,
                LinearSolverBackend::FiveBandGaussianElimination> Solver;

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

        /// @brief Fetches the current global multi-atom state vector synchronized to the Schrödinger picture.
        /// @details Evaluates the timestamp of the internal interaction picture representation. If it is newer, 
        /// the function executes an inverse transformation using the complete multi-qubit global RWA energies 
        /// to synchronize and return the physical state vector |Ψ(t)⟩.
        /// @return Const reference to the up-to-date Schrödinger state vector.
        const auto& getStateVector()
        {
            if (m_GlobalStateVectorDirac.m_TimeStamp > m_GlobalStateVectorSchrodinger.m_TimeStamp)
            {
                static const eigenenergies_t Energies = m_RwaFrame.generateGlobalRwaEnergies<QubitCount>();
                m_GlobalStateVectorSchrodinger =
                    PictureConverter::toSchrodingerPicture(m_GlobalStateVectorDirac, Energies);
            }
            return m_GlobalStateVectorSchrodinger;
        }

        /// @brief Exposes the structural manifold properties containing hardware configurations.
        /// @return Const reference to the underlying NeutralAtomManifold descriptor.
        const auto& getManifold() const
        {
            return m_Manifold;
        }
	};
}