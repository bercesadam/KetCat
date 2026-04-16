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
    /// Dummy placeholder type to represent coordinate types in non-spatial Hilbert spaces (e.g. finite-dimensional qudit spaces).
    struct _{};

	///@brief Struct representing a finite-dimensional Hilbert space
	template<natural_t Dimension>
	struct FiniteHilbertSpace
	{
		// Dimension of the Hilbert space
		static constexpr natural_t Dim = Dimension;

        using CoordinateType = _; 
	};
	
    /// @brief Enumeration to denote the grid spacing type in infinite, spatial Hilert spaces
    enum class GridType
    { 
        Uniform,
        Logarithmic
    };

	///@brief Struct representing an infinite-dimensional Hilbert space
    template<DimensionTag _SpatialDimensions, natural_t _DiscretizationSteps, real_t _SystemExtent, GridType _Grid = GridType::Uniform>
    struct InfiniteHilbertSpace
    {
        /// Expose template parameters
        static constexpr natural_t SpatialDimensions = _SpatialDimensions.value;
        static constexpr natural_t Steps = _DiscretizationSteps;
        static constexpr real_t Extent = _SystemExtent;
        static constexpr GridType Grid = _Grid;

        // Physical coordinate of grid index 0 in case of Logarithmic grids to avoid singularity at 0
        static constexpr real_t RMin = 0.001;

        /// Define type for coordinate-based indexing
        using CoordinateType = coordinate_t<_SpatialDimensions.value>; 

        /// Total number of grid points, cf. the size of the state vectors
        static constexpr natural_t Dim = ConstexprMath::pow(Steps, SpatialDimensions);

        /// @brief Uniform grid spacing Δx = Extent / N.
        /// Only available when Grid == GridType::Uniform.
        /// @details This function replaces the static constexpr data member from the earlier versions
        /// @warning Use gridSpacing(i) for grid-type-agnostic code.
        static consteval real_t dx() noexcept
            requires (Grid == GridType::Uniform)
        {
            return static_cast<real_t>(Extent) / static_cast<real_t>(Steps);
        }

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

        /// @brief Convert a per-axis grid index to a physical radius/coordinate.
        ///
        /// Uniform:      r_i = i * dx
        /// Logarithmic:  r_i = RMin * (Extent / RMin)^(i / (N-1))
        ///
        /// @param i  Per-axis grid index in [0, Steps-1].
        /// @return   Physical coordinate >= 0.
        static constexpr real_t gridToR(natural_t i) noexcept
        {
            if constexpr (Grid == GridType::Uniform)
            {
                return static_cast<real_t>(i) * dx();
            }
            else
            {
                real_t NormalizedGridPosition  = static_cast<real_t>(i) /
                           static_cast<real_t>(Steps - 1);
                return RMin * ConstexprMath::pow(Extent / RMin, NormalizedGridPosition);
            }
        }

        /// @brief Convert a physical coordinate to the nearest per-axis grid index.
        ///
        /// Uniform:      i = floor(r / dx)
        /// Logarithmic:  i = floor( log(r / RMin) / log(Extent / RMin) * (N-1) )
        ///
        /// @param r  Physical coordinate >= 0.
        /// @return   Nearest grid index, clamped to [0, Steps-1].
        static constexpr natural_t rToGrid(real_t r) noexcept
        {
            if constexpr (Grid == GridType::Uniform)
            {
                if (r <= 0.0)    return 0;
                if (r >= Extent) return Steps - 1;
                return static_cast<natural_t>(r / dx());
            }
            else
            {
                if (r <= RMin)  return 0;
                if (r >= Extent) return Steps - 1;

                real_t NormalizedGridPosition = ConstexprMath::log(r / RMin) /
                           ConstexprMath::log(Extent / RMin);
                return static_cast<natural_t>(NormalizedGridPosition * static_cast<real_t>(Steps - 1));
            }
        }

        
        /// @brief Returns the grid spacing Δr at per-axis index i.
        ///
        /// Uniform:      Δr = dx  (constant)
        /// Logarithmic:  Δr_i = r_i * (q - 1)
        ///               where q = (Extent / RMin)^(1 / (N-1))
        ///
        /// @param i  Per-axis grid index in [0, Steps-1].
        /// @return   Local grid spacing at index i.
        static constexpr real_t gridSpacing(natural_t i) noexcept
        {
            if constexpr (Grid == GridType::Uniform)
            {
                return dx();
            }
            else
            {
                real_t LogGridRatio = ConstexprMath::pow(
                    Extent / RMin,
                    1.0 / static_cast<real_t>(Steps - 1)
                );
                return gridToR(i) * (LogGridRatio - 1.0);
            }
        }

        /// @brief Returns the discrete cell hypervolume ΔV at a flat grid index.
        ///
        /// Uniform:      ΔV = dx^D  (constant, same as original discreteVolume())
        /// Logarithmic:  ΔV = Π_d Δr_{i_d}
        ///               where i_d is the per-axis index for dimension d,
        ///               recovered by decomposing the flat index.
        ///
        /// This is the integration weight for the discrete approximation:
        ///     ∫ f(x) d^D x  ≈  Σ_i f_i · ΔV_i
        ///
        /// @param flatIndex  Flattened grid index in [0, Dim-1].
        /// @return           Cell hypervolume at that grid point.
        static constexpr real_t cellVolume(natural_t flatIndex) noexcept
        {
            if constexpr (Grid == GridType::Uniform)
            {
                // Original behavior: ΔV = dx^D
                real_t Volume = 1.0;
                for (natural_t d = 0; d < SpatialDimensions; ++d)
                {
                    Volume *= dx();
                }
                return Volume;
            }
            else
            {
                // Decompose flat index → per-axis indices,
                // then multiply per-axis spacings.
                real_t    Volume    = 1.0;
                natural_t Remaining = flatIndex;

                for (natural_t d = 0; d < SpatialDimensions; ++d)
                {
                    natural_t Index = Remaining % Steps;
                    Remaining    /= Steps;
                    Volume       *= gridSpacing(Index);
                }
                return Volume;
            }
        }

    };

	///@brief Concept matches any Hilbert space
	template <typename T>
	concept hilbert_space_t =
		requires {
			{ T::Dim } -> std::convertible_to<natural_t>;	
	};

    ///@brief Concept matches any spatial (Infinite) Hilbert space
    template <typename T>
    concept spatial_hilbert_space_t =
        requires {
            { T::SpatialDimensions } -> std::convertible_to<natural_t>;
            { T::gridSpacing(natural_t{0}) } -> std::convertible_to<real_t>;
    };

    ///@brief Concept matches a spatial (Infinite) Hilbert space
    ///       of exact spatial dimensions
    template <typename T, DimensionTag _SpatialDimensions>
    concept spatial_hilbert_space_with_dim_t =
        requires {
            { T::Dim } -> std::convertible_to<natural_t>;
            requires T::SpatialDimensions == _SpatialDimensions.value;
    };

}