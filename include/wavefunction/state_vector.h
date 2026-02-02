#pragma once
#include "core_types.h"

namespace KetCat
{
	/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
	/// @tparam HilbertDim  Dimension of the Hilbert space (number of basis states).
	template <dimension_t HilbertDim>
	struct StateVector
	{
		/// Underlying state vector array
		state_vector_t<HilbertDim> m_StateVector;

	public:
		/// @brief Indexing operator
		/// @return Reference to a complex number at the given state index
		constexpr cplx_t& operator[](dimension_t index) noexcept
		{
			return m_StateVector.at(index);
		}

		/// @brief Indexing operator (const)
		/// @return Const reference to a complex number at the given state index
		constexpr const cplx_t& operator[](dimension_t index) const noexcept
		{
			return m_StateVector.at(index);
		}

		/// @brief Get the probabilities of measuring the selected basis states.
		constexpr probability_vector_t<HilbertDim> getProbabilities() const noexcept
		{
			probability_vector_t<HilbertDim> Probabilities;

			for (int i = 0; i < HilbertDim; ++i)
			{
				Probabilities[i] = m_StateVector[i].normSquared();
			}

			return Probabilities;
		}

		/// @brief  Normalize a discrete wavefunction on a uniform spatial grid
		///         so that ∑|ψᵢ|² · Δx = 1.
		/// 
		/// @tparam Dim     Dimension (number of grid points).
		/// 
		/// @param  psi     State vector representing ψ(x) or u(r) sampled on a 1D grid.
		/// @param  dx      Grid spacing Δx (or Δr for radial problems).
		/// 
		/// @details Continuous quantum wavefunctions must satisfy the normalization condition:
		///        ∫|ψ(x)|² dx = 1.
		/// 
		/// On a uniform discretized grid xᵢ = i·Δx, this integral becomes the Riemann sum:
		///        ∑|ψᵢ|² · Δx ≈ 1.
		/// 
		/// Therefore the correct discrete normalization requires multiplying the sum of
		/// squared magnitudes by Δx before taking the square root.
		/// 
		/// This routine computes:
		///        norm² = ∑ |ψᵢ|²,
		///        norm² ← norm² · Δx,
		///        ψᵢ ← ψᵢ / √(norm²).
		/// 
		constexpr void normalize(double dx) noexcept
		{
			double Norm2 = 0.0;

			// Accumulate Σ |ψᵢ|²  
			for (const cplx_t& c : m_StateVector)
			{
				Norm2 += c.normSquared();
			}

			// Convert into discrete integral: Σ |ψᵢ|² · Δx
			Norm2 *= dx;

			// Guard against division by zero
			if (Norm2 > 0.0)
			{
				const double Inv = 1.0 / ConstexprMath::sqrt(Norm2);

				// Rescale all amplitudes so that Σ |ψᵢ|² · Δx = 1
				for (cplx_t& c : m_StateVector)
				{
					c = c * Inv;
				}
			}
		}

		/// @brief Multiply this state vector by a matrix.
		/// @param mat  The matrix to multiply with (HilbertDim x HilbertDim).
		/// @return The resulting state vector.
		constexpr StateVector<HilbertDim> matMul(const matrix_t<HilbertDim>& mat) const noexcept
		{
			StateVector<HilbertDim> Result;
			for (dimension_t i = 0; i < HilbertDim; ++i)
			{
				cplx_t Sum = cplx_t::zero();
				for (dimension_t j = 0; j < HilbertDim; ++j)
				{
					Sum = Sum + mat[i][j] * m_StateVector[j];
				}
				Result[i] = Sum;
			}
			return Result;
		}
	};
}