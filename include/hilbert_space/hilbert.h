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
    struct _{};

	///@brief Struct representing a finite-dimensional Hilbert space
	template<natural_t Dimension>
	struct FiniteHilbertSpace
	{
		// Dimension of the Hilbert space
		static constexpr natural_t Dim = Dimension;

        using CoordinateType = _; 
	};
	
	///@brief Struct representing an infinite-dimensional Hilbert space
    template<DimensionTag _SpatialDimensions, natural_t _DiscretizationSteps, real_t _SystemExtent>
    struct InfiniteHilbertSpace
    {
        /// Expose template parameters
        static constexpr natural_t SpatialDimensions = _SpatialDimensions.value;
        static constexpr natural_t Steps = _DiscretizationSteps;
        static constexpr real_t Extent = _SystemExtent;

        /// Define type for coordinate-based indexing
        using CoordinateType = coordinate_t<_SpatialDimensions.value>; 

        /// Size of eg. the state vectors
        static constexpr natural_t Dim = ConstexprMath::pow(Steps, SpatialDimensions);

        /// Grid spacing Δx = Extent / N
        static constexpr real_t dx = static_cast<real_t>(Extent) / static_cast<real_t>(Steps);

        /// @brief Compute the linearized 1D index of a multidimensional coordinate
        ///        in a uniformly discretized Hilbert space.
        /// @param c Coordinate in `SpatialDimensions`‑dimensional space.
        /// @return Flattened index corresponding to `c`, assuming row‑major ordering
        ///         with `Steps` discretization points per dimensio
        static constexpr natural_t getIndex(const CoordinateType& c) noexcept {
            natural_t Index = 0;
            natural_t Multiplier = 1;

            for (natural_t i = 0; i < SpatialDimensions; ++i) {
                Index += c[i] * Multiplier;
                Multiplier *= Steps;
            }
            return Index;
        }
    };

    template<natural_t DiscretizationSteps, real_t SystemExtent>
    using InfiniteHilbertSpace1D = InfiniteHilbertSpace<1_D, DiscretizationSteps, SystemExtent>;

    template<natural_t DiscretizationSteps, real_t SystemExtent>
    using InfiniteHilbertSpace2D = InfiniteHilbertSpace<2_D, DiscretizationSteps, SystemExtent>;

	///@brief Concept to check if a type has a static member Dim of type natural_t
	template <typename T>
	concept hilbert_space_t =
		requires {
			{ T::Dim } -> std::convertible_to<natural_t>;	
	};

    template <typename T, DimensionTag _SpatialDimensions>
    concept spatial_hilbert_space_t =
        requires {
            { T::Dim } -> std::convertible_to<natural_t>;
            requires T::SpatialDimensions == _SpatialDimensions.value;
    };
}