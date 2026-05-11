#pragma once
#include "core_types.h"

namespace KetCat
{
	// Reduced Planck's constant
	// Using normalized units where ℏ = 1 for simplicity
	// In real physical units, ℏ = 1.054571817e-34 J·s
	// but for computational purposes, as the rest of the units
	// in the simulation is not in SI
	constexpr real_t hBar = 1.0;

	/// @brief Concept to ensure that the PotentialFunctor is a callable
	///        object that takes a floating-point argument and returns a floating-point value.
	/// @tparam PotentialFunctor The type of the potential functor
	/// @tparam FloatType The floating-point type used for the input and output
	template <typename PotentialFunctor, typename FloatType>
	concept potential_functor =
		std::is_floating_point_v<FloatType> &&
		std::is_class_v<PotentialFunctor> &&
		requires(PotentialFunctor obj, FloatType val)
	{
		{ obj(val) } -> std::same_as<FloatType>;
	};


	/// @brief Represents the Hamiltonian operator in 1D discreitized space
	/// @tparam Dim The dimension of the Hamiltonian matrix
	/// @details This class constructs the Hamiltonian matrix for a quantum system
	/// 		based on the provided constants and potential function.
	///			Realizes the following equation: 
	/// 		H = - (ℏ² / 2m·Δx²) · (d²/dx²) + V(x)
	template<natural_t Dim>
	class Hamiltonian
	{
		tridiagonal_matrix_t<Dim> m_hamiltonianMatrix;

	public:
		constexpr tridiagonal_matrix_t<Dim> getMatrix() const noexcept
		{
			return m_hamiltonianMatrix;
		}
		
		/// @param m   Particle mass.
		/// @param dx  Spatial discretization step Δx.
		/// @param potential  Functor representing V(x).
		/// @details The Hamiltonian is assembled as a tridiagonal matrix where:
		///          - main diagonal: 2α + V(x)
		///          - sub/super diagonals: -α
		///          with α = ℏ² / (2m Δx²).
		template<typename PotentialFunctor>
			requires potential_functor<PotentialFunctor, real_t>
		constexpr Hamiltonian(const real_t m, const real_t dx, const PotentialFunctor& potential) noexcept
		{
			// Initialize Hamiltonian matrix with zeros
			m_hamiltonianMatrix = {};

			// α = ℏ² / (2m·Δx²)
			const real_t Alpha = hBar * hBar / (2.0 * m * dx * dx);

			// Left Dirichlet point
			m_hamiltonianMatrix[MAINDIAGONAL][0] = complex_t::fromReal(1.0);

			// Right Dirichllet point
			m_hamiltonianMatrix[MAINDIAGONAL][Dim - 1] = complex_t::fromReal(1.0);

			for (natural_t i = 1; i < Dim - 1; ++i)
			{
				// Calculating position for the Potential callable: i * Δx					
				real_t Position = i * dx;

				// Superdiagonal: represents kinetic coupling to the next site (i + 1)
				if (i + 1 < Dim)
				{
					m_hamiltonianMatrix[SUPERDIAGONAL][i] = complex_t::fromReal(-Alpha);
				}

				// Main diagonal elements: Kinetic + Potential energy
				// Kinetic part: 2α (from the central term of the second-order finite difference)
				// Potential part: V(position)
				// Total: 2α + V(position)
				m_hamiltonianMatrix[MAINDIAGONAL][i] = complex_t::fromReal(2.0 * Alpha + potential(Position));

				// Subdiagonal: represents kinetic coupling to the previous site (i - 1)
				if (i > 0)
				{
					m_hamiltonianMatrix[SUBDIAGONAL][i] = complex_t::fromReal(-Alpha);
				}
			}
		}
	};
}
