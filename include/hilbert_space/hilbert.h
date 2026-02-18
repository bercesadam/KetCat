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
    template<DimensionTag _SpatialDimensions, dimension_t _DiscretizationSteps, real_t _SystemExtent>
    class InfiniteHilbertSpace
    {
    public:
        /// Define type for coordinate-based indexing for the given spatial dimensions
        /// If spatial dimensions == 1, use simple dimension_t for convenience
        using CoordinateType = std::conditional_t<
            _SpatialDimensions == 1_D,
            dimension_t,
            coordinate_t<_SpatialDimensions.value>
        >;

    private:
        static constexpr dimension_t getIndexImpl(const CoordinateType& c, dimension_t index = 0) noexcept {
            if (index == SpatialDimensions - 1)
            {
                return c[index];
            }
            else
            {
                return c[index] + _DiscretizationSteps * getIndex(c, ++index);
            }
        }

        static constexpr dimension_t calculateContainerSize() noexcept
        {
            dimension_t Size = _SpatialDimensions.value;
            for (dimension_t i = 0; i < Size - 1; ++i)
            {
                Size *= _DiscretizationSteps;
            }
            return Size;
        }

    public:
        /// Expose template parameters
        static constexpr dimension_t SpatialDimensions = _SpatialDimensions.value;
        static constexpr dimension_t Steps = _DiscretizationSteps;
        static constexpr real_t Extent = _SystemExtent;

        /// Size of eg. the state vectors
        static constexpr dimension_t Dim = calculateContainerSize();

        /// Grid spacing Δx = Extent / N
        static constexpr real_t dx = _SystemExtent / static_cast<real_t>(Steps);

        static constexpr dimension_t getIndex(const CoordinateType& c) noexcept
        {
            if (SpatialDimensions == 1)
            {
                return c;
            }
            return getIndexImpl(c);
        }
    };

    template<dimension_t DiscretizationSteps, real_t SystemExtent>
    using InfiniteHilbertSpace1D = InfiniteHilbertSpace<1_D, DiscretizationSteps, SystemExtent>;

    template<dimension_t DiscretizationSteps, real_t SystemExtent>
    using InfiniteHilbertSpace2D = InfiniteHilbertSpace<2_D, DiscretizationSteps, SystemExtent>;

	///@brief Concept to check if a type has a static member Dim of type dimension_t
	template <typename T>
	concept hilbert_space_t =
		requires {
			{ T::Dim } -> std::convertible_to<dimension_t>;	
	};

    template <typename T, DimensionTag _SpatialDimensions>
    concept spatial_hilbert_space_t =
        requires {
            { T::Dim } -> std::convertible_to<dimension_t>;
            requires T::SpatialDimensions == _SpatialDimensions.value;
    };
}