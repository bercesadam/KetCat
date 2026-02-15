#pragma once
#include <type_traits>

namespace KetCat
{
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
    template <hilbert_space_t Space>
    constexpr void normalize1D(StateVector<Space>& psi) noexcept
    {
        double Norm2 = 0.0;

        // Accumulate Σ |ψᵢ|²  
        for (const cplx_t& c : psi)
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
            for (cplx_t& c : psi)
            {
                c = c * Inv;
            }
        }
    }

    /// @brief Normalize a discrete 2D polar wavefunction on a uniform (r, θ) grid
    ///        so that ∑ |ψᵢ|² · rᵢ · Δr · Δθ = 1.
    /// 
    /// @tparam RadialSteps   Number of radial steps (from Hilbert space)
    /// @tparam AngularSteps  Number of angular steps (from Hilbert space)
    /// 
    /// @param sv             State vector representing ψ(r, θ)
    template <hilbert_space_t Space>
    constexpr void normalizePolar2D(StateVector<Space>& psi) noexcept
    {
        real_t Norm2 = 0.0;
        
        using Space = std::remove_reference_t<decltype(psi)>::HilbertSpaceType;

        // Loop over all radial and angular indices
        for (dimension_t r_idx = 0; r_idx < Space::RadDim; ++r_idx)
        {
            for (dimension_t theta_idx = 0; theta_idx < Space::AngDim; ++theta_idx)
            {
                // Convert 2D indices to 1D
                const dimension_t idx = Space::index(r_idx, theta_idx);

                // Radial coordinate (center of the cell)
                const double r = (r_idx + 0.5) * Space::dr;

                // Accumulate |ψ|² * r
                Norm2 += psi[idx].normSquared() * r;
            }
        }

        // Multiply by dr * dθ to approximate ∫ |ψ|² r dr dθ
        Norm2 *= (Space::dr * Space::dTheta);

        // Guard against division by zero
        if (Norm2 > 0.0)
        {
            const double Inv = 1.0 / ConstexprMath::sqrt(Norm2);

            for (auto& c : psi)
            {
                c = c * Inv;
            }
        }
    }
}