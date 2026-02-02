#pragma once
#include "core_types.h"
#include "constexprmath/constexpr_trigon.h"

// Aggressively include all relevant headers for ease of use
#include "hamiltonian/hamiltonian.h"
#include "hamiltonian/potential_barrier.h"
#include "hamiltonian/softculomb_potential.h"
#include "hamiltonian/harmonic_osc_potential.h"

#include "wavefunction/gaussian_wave_packet.h"
#include "wavefunction/hydrogen.h"
#include "wavefunction/eigenstate.h"
#include "wavefunction/coherent_state_gaussian.h"

#include "solvers/crank_nicolson_solver.h"


namespace KetCat
{
	/// @brief Configuration parameters for a one-dimensional particle in a box system.
	/// @tparam SpatialDiscretizationStep  Number of spatial discretization steps (including boundaries).
	/// @details As a design decision, this configuration struct was separated from the main quantum system class
	/// as the classes which are used to construct an 1D box system (initial wavefunction, Hamiltonian) often require knowledge of
	/// these parameters. This ensures that the same values are used consistently across all classes.
	template<dimension_t SpatialDiscretizationStep>
	struct OneDimensionalParticleBoxConfig
	{
		// Number of spatial discretization steps
		const dimension_t N = SpatialDiscretizationStep;
		// Dirichlet boundary conditions
		const dimension_t M = N - 2;

		// Box length in meters
		float_t L;
		// Time step in seconds
		float_t dt;
		// Spatial discretization step in meters
		float_t dx;

		/// @brief Constructs a configuration for a one-dimensional particle in a box system.
		/// @param boxLength    Length of the box.
		/// @param timeStep     Time step size for evolution.
		constexpr OneDimensionalParticleBoxConfig(float_t boxLength, float_t timeStep)
			: L(boxLength), dt(timeStep), dx(boxLength / (SpatialDiscretizationStep - 1))
		{
		}
	};

	/// @brief One-dimensional particle in a box quantum system.
	/// @tparam SpatialDiscretizationStep  Number of spatial discretization steps (including boundaries).
	template<dimension_t SpatialDiscretizationStep>
	class OneDimensionalParticleBox
	{
		//@brief Dimension of the state vector excluding boundary points
		static constexpr dimension_t StateVectorSize = SpatialDiscretizationStep - 2;

		//@brief Configuration of the particle in a box system
		OneDimensionalParticleBoxConfig<SpatialDiscretizationStep> m_config;

		//@brief Hamiltonian operator of the system
		Hamiltonian<StateVectorSize> m_hamiltonian;

		//@brief State vector of the system containing only the inner points (Dirichlet BCs)
		StateVector<InfiniteHilbertSpace<StateVectorSize>> m_psi;

		//@brief Crank-Nicolson solver for time evolution
		CrankNicolsonSolver<StateVectorSize> m_timeEvolutionSolver;

	public:
		/// @brief Constructs a one-dimensional particle in a box system.
		/// @param config        Configuration parameters for the system.
		/// @param hamiltonian   Hamiltonian operator of the system.
		/// @param stateVector   Initial state vector of the system.
		constexpr OneDimensionalParticleBox(
			const OneDimensionalParticleBoxConfig<SpatialDiscretizationStep>& config,
			const Hamiltonian<StateVectorSize>& hamiltonian, const StateVector<InfiniteHilbertSpace<StateVectorSize>>& stateVector) noexcept
			: m_config(config), m_hamiltonian(hamiltonian), m_psi(stateVector),
			m_timeEvolutionSolver(hamiltonian, config.dt)
		{
		}

		/// @brief Evolves the system by one time step using the Crank-Nicolson method.
		constexpr StateVector<InfiniteHilbertSpace<StateVectorSize>> evolve() noexcept
		{
			m_psi = m_timeEvolutionSolver(m_psi);
			return m_psi;
		}
	};
} // namespace KetCat
	