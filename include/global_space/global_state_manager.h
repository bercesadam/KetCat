#pragma once
#include <bitset>
#include <algorithm> 
#include <ranges>

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

		/// @brief Full system state vector in the product Hilbert space.
		StateVector<typename SubspaceManager::FullHilbertSpace> m_GlobalStateVector{};

    public:
		GlobalStateManager(std::bitset<QubitCount> initialState)
		{
			m_GlobalStateVector =
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

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(HartreeEnergies, DipoleMatrix, lasers);

            static CrankNicolsonSolver<ConfigType::LevelCount, LinearSolverBackend::ThomasTridiagonal> Solver;

            Hamiltonian.updateMainDiagonal(lasers);
            Hamiltonian.updateOffDiagonal(lasers);
            Solver.updateMatrices(Hamiltonian.getMatrix(), TimeMaster::Clock().getTimeStep());

            // Map the local 1-qubit Hamiltonian operation to the global N-qubit state vector
            std::array<natural_t, 1> targets = { targetAtom };
            SubspaceManager::template performTimeEvolution<1>(Solver, m_GlobalStateVector, targets);
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

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                SingleAtomExcitation(HartreeEnergies, DipoleMatrix, lasers);

            static TwoAtomRydbergBlockade<ConfigType::LevelCount>
                RydbergBlockade(Units::MeterToAtomicLength * 1E-9,
                    ConfigType::RydbergLevel, HartreeEnergies, DipoleMatrix);

            static CrankNicolsonSolver<ConfigType::LevelCount,
                LinearSolverBackend::FiveBandGaussianElimination> Solver;

            SingleAtomExcitation.updateMainDiagonal(lasers);
            SingleAtomExcitation.updateOffDiagonal(lasers);
            auto SingleAtomHamiltonian = SingleAtomExcitation.getMatrix();

            RydbergBlockade.updateMatrix(SingleAtomHamiltonian, SingleAtomHamiltonian);

            // --- DIRAC (INTERACTION) PICTURE TRANSFORMATION ---
            auto H_Schrodinger = RydbergBlockade.getMatrix();
            auto H_Interaction =
                InteractionPictureHamiltonian<ConfigType::LevelCount>::transform
                (H_Schrodinger, TimeMaster::Clock().getGlobalTime());

            Solver.updateMatrices(H_Interaction, TimeMaster::Clock().getTimeStep());

            std::array<natural_t, 2> targets = { controlAtom, targetAtom };
            SubspaceManager::template performTimeEvolution<2>(Solver, m_GlobalStateVector, targets);
        }

        const auto& getStateVector() const
        {
            return m_GlobalStateVector;
        }

        const auto& getManifold() const
        {
            return m_Manifold;
        }
	};
}