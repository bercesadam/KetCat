#pragma once

#include "core_types.h"


namespace KetCat
{
     /// @brief 2D Cartesian coordinate Hilbert space.
     /// @details
     ///   Uniform discretization on a square domain.
     ///   The physical domain is assumed to span
     ///       [-Extent/2, +Extent/2]
     ///   along both axes.
     /// 
     /// @tparam DiscretizationSteps Number of grid points per axis.
     /// @tparam SystemExtent Total physical length of the domain.
    template<dimension_t DiscretizationSteps, real_t SystemExtent>
    struct Hilbert2D
    {
        static constexpr dimension_t Steps = DiscretizationSteps;

        /// Total Hilbert space dimension.
        static constexpr dimension_t Dim =
            DiscretizationSteps * DiscretizationSteps;

        /// Grid spacing Δx = Extent / N
        static constexpr real_t dx =
            SystemExtent / DiscretizationSteps;

        /// Converts (ix, iy) → linear index.
        static constexpr dimension_t
            getIndex(dimension_t ix, dimension_t iy) noexcept
        {
            return ix + DiscretizationSteps * iy;
        }
    };
}
