#pragma once
#include "hilbert_space/state_vector.h"
#include "hilbert_space/hilbert.h"
#include "constexprmath/constexpr_core_functions.h"

namespace KetCat
{
    /// @brief Compute physicists’ Hermite polynomial Hₙ(x) using recursion.
    /// @details
    /// Recurrence:
    ///   H₀(x) = 1  
    ///   H₁(x) = 2x  
    ///   Hₙ(x) = 2x·Hₙ₋₁(x) − 2(n−1)·Hₙ₋₂(x)  for n ≥ 2
    constexpr real_t hermite(natural_t n, real_t x) noexcept
    {
        if (n == 0) return 1.0;
        if (n == 1) return 2.0 * x;

        real_t H0 = 1.0;
        real_t H1 = 2.0 * x;
        real_t Hn = 0.0;

        for (natural_t i = 2; i <= n; ++i)
        {
            Hn = 2.0 * x * H1 - 2.0 * (i - 1) * H0;
            H0 = H1;
            H1 = Hn;
        }
        return Hn;
    }

    /// @brief 2D quantum harmonic oscillator seed wavefunction generator.
    /// @details
    /// Constructs a Hermite–Gaussian seed of the form:
    ///
    ///     ψₙₓ,ₙᵧ(x,y) = Hₙₓ(αx) · Hₙᵧ(αy) · exp(−½ α² (x² + y²))
    ///
    /// where:
    ///   • Hₙ  = physicists’ Hermite polynomial  
    ///   • α   = √(m·ω) scaling factor (atomic units)
    ///
    /// In atomic units the harmonic oscillator width is:
    ///
    ///     a₀ = 1 / α
    ///
    /// For trapped ions (e.g. ⁴⁰Ca⁺), α ≪ 1 because the mass is large.
    ///
    /// When α corresponds to the real ion mass and trap frequency,
    /// the generated wavefunctions match physical motional states used
    /// in trapped‑ion quantum computing.
    template<spatial_hilbert_space_with_dim_t<2_D> HilbertSpace>
    struct Harmonic2D
    {
        /// @brief Generate ψₙₓ,ₙᵧ(x,y) on the Hilbert-space grid.
        /// @details
        /// The returned wavefunction is *not* normalized.
        /// Use StateVector.normalize() to enforce:
        ///
        ///     Σ |ψᵢⱼ|² (Δx)² = 1.
        ///
        /// @param nx  Quantum number nₓ (0,1,2,…)
        /// @param ny  Quantum number nᵧ (0,1,2,…)
        /// @param alpha Scaling parameter α controlling wavefunction width
        ///              α = √(m·ω)
        ///               • m: particle mass (a.u.)
        ///               • ω: angular trap frequency (a.u.)
        ///              Examples:
        ///               • Electron: α ≈ 1  
        ///           • Ca⁺ in 1 MHz trap: α ≈ 0.003346  
        /// @return    A StateVector with ψₙₓ,ₙᵧ sampled on the grid.
        StateVector<HilbertSpace> operator()(natural_t nx, natural_t ny, real_t alpha)
        {
            StateVector<HilbertSpace> Psi{ cplx_t::zero() };

            // Grid spacing Δx = Extent / (Dim-1)
            const real_t dx = HilbertSpace::Extent / (HilbertSpace::Dim - 1);
            constexpr natural_t Steps = HilbertSpace::Steps;

            for (natural_t ix = 0; ix < Steps; ++ix)
            {
                for (natural_t iy = 0; iy < Steps; ++iy)
                {
                    // --- 1) Convert grid indices → continuous coordinates (centered)
                    real_t x = (real_t(ix) - Steps * 0.5) * dx;
                    real_t y = (real_t(iy) - Steps * 0.5) * dx;

                    // --- 2) Evaluate Hermite polynomials Hₙ(α x)
                    real_t Hx = hermite(nx, alpha * x);
                    real_t Hy = hermite(ny, alpha * y);

                    // --- 3) Gaussian envelope exp(−½ α² (x²+y²))
                    real_t Envelope = 
                        ConstexprMath::exp<20>(-0.5 * alpha * alpha * (x*x + y*y));

                    // --- 4) Full 2D Hermite–Gauss mode
                    cplx_t Value = cplx_t(Hx * Hy * Envelope, 0.0);

                    Psi[{ ix, iy }] = Value;
                }
            }
            
            Psi.normalize();
            return Psi;
        }
    };

}