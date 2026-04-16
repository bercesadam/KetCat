#pragma once
#include "hilbert.h"

namespace KetCat
{
	/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
	/// @tparam Dim  Dimension of the Hilbert space (number of basis states).
	template <hilbert_space_t HilbertSpace>
	struct StateVector
	{
		/// Type alias for the Hilbert space type
		using HilbertSpaceType = HilbertSpace;

		/// Size of the vector for convenience
		static constexpr natural_t Size = HilbertSpace::Dim;

		/// Underlying state vector array
		state_vector_t<Size> m_StateVector;

		/// INDEXING ///////////////////////////////////////////////////////////////////////

		/// @brief Indexing operator using the coordinate type (array of natural_t's)
		///        provided by the Hilbert Space type (depending on the spatial dimensions)
		/// @param Spatial coordinates on the discrete grid
		/// @return Reference to a complex number at the given state index
		constexpr complex_t& operator[](HilbertSpace::CoordinateType index) noexcept
		{
			return m_StateVector.at(HilbertSpace::getIndex(index));
		}

		/// @brief Const indexing operator using the coordinate type (array of natural_t's)
		///       provided by the Hilbert Space type (depending on the spatial dimensions)
		/// @param Spatial coordinates on the discrete grid
		/// @return Const reference to a complex number at the given state index
		constexpr const complex_t& operator[](HilbertSpace::CoordinateType index) const noexcept
		{
			return m_StateVector.at(HilbertSpace::getIndex(index));
		}

		/// @brief Indexing operator to access the flattened state vector,
		///        regardless the spatial dimensions of the underlying Hilbert space.
		/// @param Index of the discrete grid
		/// @return Reference to a complex number at the given state index
		constexpr complex_t& operator[](natural_t index) noexcept
		{
			return m_StateVector.at(index);
		}

		/// @brief Const indexing operator to access the flattened state vector,
		///        regardless the spatial dimensions of the underlying Hilbert space.
		/// @param Index of the discrete grid
		/// @return Const reference to a complex number at the given state index
		constexpr const complex_t& operator[](natural_t index) const noexcept
		{
			return m_StateVector.at(index);
		}

		/// BASIC LINEAR ALGEBRA ///////////////////////////////////////////////////////////////

		/// @brief Multiply this state vector by a matrix.
		/// @param mat  The matrix to multiply with (Size x Size).
		/// @return The resulting state vector.
		constexpr StateVector<HilbertSpaceType> matMul(const matrix_t<Size>& mat) const noexcept
		{
			StateVector<HilbertSpaceType> Result;
			for (natural_t i = 0; i < Size; ++i)
			{
				complex_t Sum = complex_t::zero();
				for (natural_t j = 0; j < Size; ++j)
				{
					Sum = Sum + mat[i][j] * m_StateVector[j];
				}
				Result[i] = Sum;
			}
			return Result;
		}

		/// @brief Compute the inner product ⟨ψ|φ⟩.
		///
		/// Continuous definition:
		///     ⟨ψ|φ⟩ = ∫ ψ*(x) φ(x) d^D x
		///
		/// Discrete approximation:
		///     ⟨ψ|φ⟩ ≈ Σ_i ψ_i* φ_i · ΔV
		///
		/// where:
		///     ψ_i*  = complex conjugate of ψ_i
		///     ΔV    = dx^D
		///
		/// @param phi    The state |φ⟩.
		/// @return       The complex amplitude ⟨ψ|φ⟩.
		constexpr complex_t innerProduct(const StateVector<HilbertSpaceType>& phi) const noexcept
			requires (spatial_hilbert_space_t<HilbertSpaceType>)
		{
			complex_t Result = complex_t::zero();

			for (natural_t i = 0; i < Size; ++i)
			{
				// ψ_i* · φ_i
				Result += m_StateVector[i].conj() * phi.m_StateVector[i] * HilbertSpaceType::cellVolume(i);;
			}

			// Multiply by discrete cell volume ΔV = dx^D
			return Result;
		}

		/// PROBABILITY MEASUREMENTS ////////////////////////////////////////////////////////

		/// @brief Compute the probability of projecting onto another state |φ⟩.
		///
		/// Born rule:
		///     P = |⟨ψ|φ⟩|²
		///
		/// Using the discrete inner product:
		///     ⟨ψ|φ⟩ = Σ_i ψ_i* φ_i · ΔV
		///
		/// @param phi    The state |φ⟩.
		/// @return       Probability in range [0, 1].
		constexpr real_t probabilityOf(const StateVector<HilbertSpaceType>& phi) const noexcept
			requires (spatial_hilbert_space_t<HilbertSpaceType>)
		{
			const complex_t Amplitude = innerProduct(phi);
			return Amplitude.normSquared();
		}

		/// @brief Get the probabilities of measuring the selected basis states.
		/// @returns Array of real_t's with the probability of each spatial position
		constexpr probability_vector_t<Size> getProbabilities() const noexcept
		{
			probability_vector_t<Size> Probabilities;

			for (natural_t i = 0; i < Size; ++i)
			{
				Probabilities[i] = m_StateVector[i].normSquared();
			}

			return Probabilities;
		}

		/// SPATIAL DIMENSION-AGNOSTIC NORMALIZATION ////////////////////////////////////////
		
		/// @brief Normalize the state vector so that the total probability equals 1.
		///
		/// Continuous normalization condition:
		///     ∫ |ψ(x)|² d^D x = 1
		///
		/// Discrete approximation on a uniform grid:
		///     Σ_i |ψ_i|² · ΔV = 1
		///
		/// where:
		///     ΔV = dx^D  (cell hypervolume)
		///
		/// Procedure:
		///     1) Compute norm² = Σ_i |ψ_i|² · ΔV
		///     2) Rescale:
		///            ψ_i ← ψ_i / √(norm²)
		///
		/// After normalization:
		///     ⟨ψ|ψ⟩ = 1
		constexpr void normalize() noexcept
			requires (spatial_hilbert_space_t<HilbertSpaceType>)
		{
			// Compute ⟨ψ|ψ⟩ = Σ ψ_i* ψ_i · ΔV
			const real_t Norm2 = innerProduct(*this).re;

			if (Norm2 > 0.0)
			{
				const real_t InvNorm = 1.0 / ConstexprMath::sqrt(Norm2);

				for (complex_t& c : m_StateVector)
					c = c * InvNorm;
			}
		}

		/// MANUAL SUPERPOSITION ////////////////////////////////////////////////////////////

		/// @brief Create a superposition of two state vectors with given coefficients.
		/// @param phi    The other state vector to superpose with.
		/// @param alpha  Complex amplitude for this state vector.
		/// @param beta   Complex amplitude for the other state vector.
		/// @param dx     Grid spacing for normalization.
		constexpr StateVector<HilbertSpaceType> superpose(const StateVector<HilbertSpaceType>& phi,
			complex_t alpha, complex_t beta) const noexcept
			requires (spatial_hilbert_space_t<HilbertSpaceType>)
		{
			StateVector Result;

			for (natural_t i = 0; i < StateVector::Size; ++i)
			{
				// |Ψ⟩ = α |ψ₀⟩ + β |ψ₁⟩
				Result[i] = alpha * m_StateVector[i] + beta * phi.m_StateVector[i];
			}

			Result.normalize();
			return Result;
		}
	};
}
