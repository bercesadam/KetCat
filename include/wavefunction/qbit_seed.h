#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"

namespace KetCat::QCC
{
	/// @brief Functor to generate computational basis states for QBitCount qubits.
	/// @tparam QBitCount  Number of qubits.
	template<natural_t QBitCount>
	struct QBitState
	{
		/// @brief Number of basis states (2^QBitCount)
		static constexpr natural_t Dim = ConstexprMath::pow2(QBitCount);

		/// @brief Generates the computational basis state vector |initialValue⟩.
		/// @param initialValue  Index of the basis state to generate (0 ≤ initialValue < Dim).
		/// @return State vector representing the selected basis state.
		constexpr StateVector<FiniteHilbertSpace<Dim>>
			operator()(unsigned int initialValue = 0) const noexcept
		{
			StateVector<FiniteHilbertSpace<Dim>> Psi{};
			Psi[initialValue] = cplx_t::fromReal(1.0);
			return Psi;
		}
	};
}