#pragma once
#include "hilbert.h"

namespace KetCat
{
	/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
	/// @tparam Dim  Dimension of the Hilbert space (number of basis states).
	template <hilbert_space_t Space>
	struct StateVector
	{
		/// Size of the vector for convenience
		static constexpr dimension_t Size = Space::Dim;


		/// Type alias for the Hilbert space type
		using HilbertSpaceType = Space;

		/// Underlying state vector array
		state_vector_t<Size> m_StateVector;


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
		constexpr probability_vector_t<Size> getProbabilities() const noexcept
		{
			probability_vector_t<Size> Probabilities;

			for (int i = 0; i < Size; ++i)
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
		constexpr void normalize(real_t dx) noexcept
		{
			real_t Norm2 = 0.0;

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
				const real_t Inv = 1.0 / ConstexprMath::sqrt(Norm2);

				// Rescale all amplitudes so that Σ |ψᵢ|² · Δx = 1
				for (cplx_t& c : m_StateVector)
				{
					c = c * Inv;
				
				}
			}
		}

		constexpr void normalize2D(real_t dx) noexcept
		{
			real_t Norm2 = 0.0;

			for (const cplx_t& c : m_StateVector)
				Norm2 += c.normSquared();

			// Diszkrét 2D integrál: Σ |ψ|² · dx²
			Norm2 *= Space::dx * dx;

			if (Norm2 > 0.0)
			{
				const real_t Inv = 1.0 / ConstexprMath::sqrt(Norm2);
				for (cplx_t& c : m_StateVector)
					c = c * Inv;
			}
		}

		/// @brief Multiply this state vector by a matrix.
		/// @param mat  The matrix to multiply with (Size x Size).
		/// @return The resulting state vector.
		constexpr StateVector<HilbertSpaceType> matMul(const matrix_t<Size>& mat) const noexcept
		{
			StateVector<HilbertSpaceType> Result;
			for (dimension_t i = 0; i < Size; ++i)
			{
				cplx_t Sum = cplx_t::zero();
				for (dimension_t j = 0; j < Size; ++j)
				{
					Sum = Sum + mat[i][j] * m_StateVector[j];
				}
				Result[i] = Sum;
			}
			return Result;
		}

		/// @brief Create a superposition of two state vectors with given coefficients.
		/// @param other  The other state vector to superpose with.
		/// @param alpha  Complex amplitude for this state vector.
		/// @param beta   Complex amplitude for the other state vector.
		/// @param dx     Grid spacing for normalization.
		constexpr StateVector<HilbertSpaceType> superpose(const StateVector<HilbertSpaceType>& other,
			cplx_t alpha, cplx_t beta) const noexcept
		{
			StateVector Result;

			for (dimension_t i = 0; i < StateVector::Size; ++i)
			{
				// |Ψ⟩ = α |ψ₀⟩ + β |ψ₁⟩
				Result[i] = alpha * m_StateVector[i] + beta * other.m_StateVector[i];
			}

			//Result.normalize(dx);

			return Result;
		}

		/// @brief Compute the inner product ⟨ψ|φ⟩ between this and another state vector.
		/// @param other  The other state vector |φ⟩.
		/// @param dx     Grid spacing for proper normalization.
		/// @return The inner product ⟨ψ|φ⟩ as a complex number.
		constexpr cplx_t innerProduct(const StateVector<HilbertSpaceType>& other, double dx) const noexcept
		{
			cplx_t Result{ cplx_t::zero() };

			for (dimension_t i = 0; i < StateVector::Size; ++i)
			{
				// ⟨ψ|φ⟩ = Σ ψᵢ* · φᵢ
				Result += m_StateVector[i].conj() * other.m_StateVector[i];
			}

			return Result * dx;
		}

		/// @brief Compute the probability of measuring the system in the state represented by another state vector.
		/// @param overlap  The other state vector |φ⟩.
		/// @param dx       Grid spacing for proper normalization.
		/// @return The probability as a float.
		constexpr real_t probabilityOf(const StateVector<HilbertSpaceType> overlap, double dx) const noexcept
		{
			const cplx_t amplitude = innerProduct(overlap, dx);
			return amplitude.normSquared();
		}
	};
}