#pragma once

#include "core_types.h"

///@file This file defines structures and concepts related to Hilbert spaces.
///      It was introduced to clearly separate the underlying concept behind
///		 the meaning of the size of the StateVector type which is used as the
/// 	 basic type in the whole KetCat library:
/// 		- Finite-dimensional Hilbert spaces for quantum computing systems.
///			- Infinite-dimensional Hilbert spaces for wavefunctions, where the
/// 			dimension represents a discretization steps of a continuous space.
///		 ...thus preventing confusion between these two different meanings of "dimension"
/// 	 and also preventing possible misuses of StateVector in the two separate contexts,
///      while allowing that both concepts share the same underlying mathematical foundation.

namespace KetCat
{
	///@brief Struct representing a finite-dimensional Hilbert space
	template<dimension_t Dimension>
	struct FiniteHilbertSpace
	{
		// Dimension of the Hilbert space
		static constexpr dimension_t Dim = Dimension;
	};
	
	///@brief Struct representing an infinite-dimensional Hilbert space
	template<dimension_t DiscretizationSteps>
	struct InfiniteHilbertSpace
	{
		// Dimension of the Hilbert space
		static constexpr dimension_t Dim = DiscretizationSteps;
	};

	///@brief Concept to check if a type has a static member Dim of type dimension_t
	template <typename T>
	concept hilbert_space_t =
		requires {
			{ T::Dim } -> std::convertible_to<dimension_t>;
	};
}