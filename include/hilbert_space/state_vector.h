#pragma once
#include "hilbert.h"

namespace KetCat
{
	/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
	/// @tparam Dim  Dimension of the Hilbert space (number of basis states).
	template <hilbert_space_t HilbertSpace>
	struct StateVector
	{
		/// Size of the vector for convenience
		static constexpr dimension_t Size = HilbertSpace::Dim;


		/// Type alias for the Hilbert space type
		using HilbertSpaceType = HilbertSpace;

		/// Underlying state vector array
		state_vector_t<Size> m_StateVector;


		/// @brief Indexing operator
		/// @return Reference to a complex number at the given state index
		constexpr cplx_t& operator[](HilbertSpace::CoordinateType index) noexcept
		{
			return m_StateVector.at(HilbertSpace::getIndex(index));
		}

		/// @brief Indexing operator (const)
		/// @return Const reference to a complex number at the given state index
		constexpr const cplx_t& operator[](HilbertSpace::CoordinateType index) const noexcept
		{
			return m_StateVector.at(HilbertSpace::getIndex(index));
		}

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

			for (dimension_t i = 0; i < Size; ++i)
			{
				Probabilities[i] = m_StateVector[i].normSquared();
			}

			return Probabilities;
		}

		/// @brief Normalize the discrete wavefunction in any spatial dimension D
		///        so that the discrete integral over the grid equals 1:
		///        
		///        Continuous normalization: ∫ |ψ(x)|² d^D x = 1
		///        Discrete approximation on a uniform grid:
		///           Σ_i |ψ_i|² · (dx)^D ≈ 1
		constexpr void normalize() noexcept
		{
			real_t Norm2 = 0.0;

			// Step 1: Accumulate the sum of squared magnitudes |ψ_i|²
			for (const cplx_t& c : m_StateVector)
				Norm2 += c.normSquared();

			// Step 2: Convert sum into a discrete D-dimensional integral: Σ |ψ_i|² · dx^D
			constexpr dimension_t D = HilbertSpaceType::SpatialDimensions;  // Spatial dimension
			real_t dxPowD = 1.0;
			for (dimension_t d = 0; d < D; ++d)
				dxPowD *= HilbertSpaceType::dx;

			Norm2 *= dxPowD;

			// Step 3: Guard against division by zero and rescale
			if (Norm2 > 0.0)
			{
				const real_t Inv = 1.0 / ConstexprMath::sqrt(Norm2);

				// ψ_i ← ψ_i / √(Σ |ψ_i|² · dx^D)
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
