#pragma once
#include "hilbert_space/state_vector.h"
#include "hilbert_space/hilbert.h"

namespace KetCat
{
    /// @brief Compute Hermite polynomial Hₙ(x) using recursion.
    /// @details
    /// Recurrence:
    ///   H₀(x) = 1  
    ///   H₁(x) = 2x  
    ///   Hₙ(x) = 2x·Hₙ₋₁(x) − 2(n−1)·Hₙ₋₂(x)  for n ≥ 2
    ///
    /// Stable and constexpr‑friendly for moderate n.
    constexpr real_t hermite(natural_t n, real_t x) noexcept
    {
        if (n == 0)
        {
            return 1.0;
        }
        if (n == 1)
        {
            return 2.0 * x;
        }

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

    /// @brief 2D harmonic oscillator seed wavefunction generator.
    /// @details
    /// Constructs a separable product state of the form  
    ///   ψₙₓ,ₙᵧ(x,y) ∝ Hₙₓ(αx) · Hₙᵧ(αy) · exp(−½·α²·(x² + y²))  
    /// where Hₙ are physicists’ Hermite polynomials and α is a scale factor
    /// controlling the oscillator length.  
    ///
    /// This generator uses α = 1 to produce simple, shape‑correct (but not fully
    /// normalized) basis seeds for (nₓ,nᵧ) energy eigenstates of a 2D harmonic well.
    ///
    /// ---
    /// Usage in Trapped‑Ion Qubit Emulation
    ///
    /// Although this is a general mathematical 2D harmonic oscillator, its structure
    /// mirrors the motional degrees of freedom of trapped‑ion qubits:
    ///
    ///  • A single trapped ion in a radiofrequency (Paul) trap experiences an **approximately
    ///    harmonic confinement potential** in each spatial direction.  
    ///  
    ///  • The corresponding vibrational states |n⟩ of the ion’s motion are precisely the
    ///    **Hermite–Gaussian eigenfunctions** of a quantum harmonic oscillator.  
    ///  
    ///  • These vibrational levels form the “phonon modes” used in trapped‑ion quantum
    ///    gates (e.g., Mølmer–Sørensen, Cirac–Zoller), where internal qubit states couple
    ///    through shared motional excitations.  
    ///
    /// By providing ψₙₓ,ₙᵧ(x,y) for arbitrary (nₓ,nᵧ), this class can therefore be used as
    /// a **lightweight, grid‑based approximation** of trapped‑ion motional states.
    template<natural_t Dim, real_t Extent>
    struct Harmonic2D
    {
        using HilbertSpace = InfiniteHilbertSpace2D<Dim, Extent>;

        /// @brief Generate a 2D oscillator seed ψₙₓ,ₙᵧ(x,y).
        /// @param nx Quantum number along x.
        /// @param ny Quantum number along y.
        /// @return StateVector<HilbertSpace> containing ψ over the grid.
        StateVector<HilbertSpace> operator()(natural_t nx, natural_t ny)
        {
            StateVector<HilbertSpace> Psi{ cplx_t::zero() };

            const real_t dx = Extent / (Dim - 1);
            const real_t Alpha = 1.0;  // scale factor for the seed

            for (natural_t ix = 0; ix < Dim; ++ix)
            {
                for (natural_t iy = 0; iy < Dim; ++iy)
                {
                    // Center grid at origin
                    real_t x = (static_cast<real_t>(ix) - Dim / 2) * dx;
                    real_t y = (static_cast<real_t>(iy) - Dim / 2) * dx;

                    // Hermite polynomial factors
                    real_t Hx = hermite(nx, Alpha * x);
                    real_t Hy = hermite(ny, Alpha * y);

                    // Gaussian envelope exp(−½ α²(x² + y²))
                    real_t Envelope =
                        ConstexprMath::exp<20>(-0.5 * Alpha * Alpha * (x * x + y * y));

                    cplx_t value = cplx_t(Hx * Hy * Envelope, 0.0);

                    natural_t index = HilbertSpace::getIndex({ ix, iy });
                    Psi[index] = value;
                }
            }
            return Psi;
        }
    };
}
