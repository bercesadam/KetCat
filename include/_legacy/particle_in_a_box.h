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
	/// @brief One-dimensional particle in a box quantum system.
	/// @tparam SpatialDiscretizationStep  Number of spatial discretization steps (including boundaries).
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
	class OneDimensionalParticleBox
	{
		//@brief State vector of the system containing
		StateVector<HilbertSpace> m_psi;

		//@brief Crank-Nicolson solver for time evolution
		CrankNicolsonSolver<HilbertSpace> m_timeEvolutionSolver;

	public:
		/// @brief Constructs a one-dimensional particle in a box system.
		/// @param config        Configuration parameters for the system.
		/// @param hamiltonian   Hamiltonian operator of the system.
		/// @param stateVector   Initial state vector of the system.
		constexpr OneDimensionalParticleBox(
			const Hamiltonian<HilbertSpace::Dim>& hamiltonian, const StateVector<HilbertSpace>& stateVector, real_t dt) noexcept
			: m_psi(stateVector),
			  m_timeEvolutionSolver(hamiltonian.getMatrix(), dt)
		{
		}

		/// @brief Evolves the system by one time step using the Crank-Nicolson method.
		constexpr StateVector<HilbertSpace> evolve() noexcept
		{
			m_psi = m_timeEvolutionSolver(m_psi);
			return m_psi;
		}
	};
} // namespace KetCat
	