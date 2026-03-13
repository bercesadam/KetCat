#pragma once
#include "wavefunction.h"
#include "quantum_number.h"
#include "hydrogen.h"

namespace KetCat
{
    /// @brief Associated Legendre polynomial Pₗᵐ(x).
    /// @details Handles negative m using the standard Condon–Shortley phase:
    ///          Pₗ⁻ᵐ(x) = (−1)ᵐ ( (l−m)! / (l+m)! ) Pₗᵐ(x)
    inline real_t legendre(int L, int M, real_t X)
    {
        // Handle negative m using reflection formula:
        if (M < 0)
        {
            real_t Sign = (M % 2 == 0) ? 1.0 : -1.0;
            return Sign *
                ConstexprMath::factorial(L - std::abs(M)) /
                ConstexprMath::factorial(L + std::abs(M)) *
                legendre(L, -M, X);
        }

        // Compute P_M^M(x)  (diagonal term)
        real_t Pmm = 1.0;
        if (M > 0)
        {
            real_t Sqrt = ConstexprMath::sqrt((1.0 - X) * (1.0 + X));
            real_t Factor = 1.0;
            for (int i = 1; i <= M; ++i)
            {
                Pmm *= -Factor * Sqrt;
                Factor += 2.0;
            }
        }

        if (L == M)
        {
            return Pmm;
        }

        // Compute P_{M+1}^M(x)
        real_t Pm1m = X * (2 * M + 1) * Pmm;
        if (L == M + 1)
        {
            return Pm1m;
        }

        // Recurrence for higher l:
        //  Pₗᵐ = ((2l−1)x P_{l−1}ᵐ − (l+m−1)P_{l−2}ᵐ ) / (l−m)
        real_t Prev = Pmm;
        real_t Curr = Pm1m;

        for (int LL = M + 2; LL <= L; ++LL)
        {
            real_t Next = ((2 * LL - 1) * X * Curr -
                           (LL + M - 1) * Prev) / (LL - M);

            Prev = Curr;
            Curr = Next;
        }
        return Curr;
    };

    /// @brief Spherical harmonic Yₗᵐ(θ,φ).
    /// @details Normalized using the standard Nₗᵐ factor:
    ///          Yₗᵐ = Nₗᵐ Pₗᵐ(cosθ) e^{i m φ}
    inline complex_t spherical_harmonic(int L, int M, real_t Theta, real_t Phi)
    {
        real_t Norm =
            std::sqrt((2.0 * L + 1.0) /
                      (4.0 * ConstexprMath::Pi) *
                      ConstexprMath::factorial(L - ConstexprMath::abs(M)) /
                      ConstexprMath::factorial(L + ConstexprMath::abs(M)));

        real_t X = ConstexprMath::cos(Theta);
        real_t P = legendre(L, ConstexprMath::abs(M), X);

        // Phase factor e^{i m φ}
        complex_t Phase = ConstexprMath::exp<20>(complex_t(0.0, M * Phi));

        if (M < 0)
        {
            // Condon–Shortley phase
            real_t Sign = (M % 2 == 0) ? 1.0 : -1.0;
            return Phase.conj() * (Sign * Norm * P);
        }

        return Phase * (Norm * P);
    }


    // ---------------------------------------------------------------------------
    // Hydrogen 2D slice generator
    // ---------------------------------------------------------------------------
    /// @brief 2D slice of a hydrogenic wavefunction Ψₙₗₘ, fixed at y ≈ 0.
    /// @details
    /// We evaluate the full 3D formula
    ///     Ψₙₗₘ(r,θ,φ) = Rₙₗ(r) · Yₗᵐ(θ,φ)
    ///
    /// on a 2D grid (x,z) with y fixed (near 0), giving a planar cross‑section
    /// of the orbital. Produces nice visual patterns:
    ///     s orbitals → spherical blobs  
    ///     p orbitals → dumbbells  
    ///     d orbitals → clover shapes  
    ///     f orbitals → flower‑like symmetries  
    ///
    /// This is a visually clean way to inspect angular structure of orbitals.
    template<spatial_hilbert_space_with_dim_t<2_D> HilbertSpace>
    struct Hydrogen2D
    {
        /// @brief Generate a 2D hydrogenic orbital for (n,l,m).
        /// @param QNumbers QuantumNumber object containing (n,l,m).
        /// @return StateVector<HilbertSpace> containing Ψ(x,z) values.
        constexpr Wavefunction<HilbertSpace> operator()(QuantumNumber QNumbers) noexcept
        {
            constexpr natural_t Steps = HilbertSpace::Steps;
            constexpr real_t dx = HilbertSpace::dx;

            StateVector<HilbertSpace> Psi{ complex_t::zero() };

            // --------------------------------------------------------
            // 1D radial component Rₙₗ(r) = uₙₗ(r) / r
            // Already normalized by HydrogenOrbital<Dim>().
            // --------------------------------------------------------
            auto RadialArray = HydrogenOrbital<InfiniteHilbertSpace<1_D, Steps, HilbertSpace::Extent>>{}(
                QNumbers,
                1.0   // effective Bohr radius
            );

            for (natural_t ix = 0; ix < Steps; ++ix)
            {
                for (natural_t iz = 0; iz < Steps; ++iz)
                {
                    // --------------------------------------------
                    // Physical grid coordinates (centered)
                    // --------------------------------------------
                    real_t XCoord = (static_cast<real_t>(ix) - Steps / 2) * dx;
                    real_t ZCoord = (static_cast<real_t>(iz) - Steps / 2) * dx;
                    real_t YCoord = 0.01;   // Slight offset to avoid φ undefined at x=0

                    // Spherical radius
                    real_t R = ConstexprMath::sqrt(
                        XCoord * XCoord +
                        YCoord * YCoord +
                        ZCoord * ZCoord
                    );

                    complex_t PsiValue = complex_t::zero();

                    if (R > 1e-12)
                    {
                        // --------------------------------------------
                        // Radial part Rₙₗ(r) = uₙₗ(r) / r
                        // --------------------------------------------
                        natural_t IR = static_cast<natural_t>(R / dx);
                        if (IR >= Steps)
                        {
                            IR = Steps - 1;
                        }

                        real_t U = RadialArray[IR].re;
                        real_t RadialPart = U / R;

                        // --------------------------------------------
                        // Angular coordinates
                        // θ = arccos(z/r)
                        // φ = atan2(y,x)
                        // --------------------------------------------
                        real_t Theta = ConstexprMath::acos(ZCoord / R);
                        real_t Phi   = ConstexprMath::atan2(YCoord, XCoord);

                        // Angular spherical harmonic Yₗᵐ(θ,φ)
                        complex_t Ylm = spherical_harmonic(
                            QNumbers.l(),
                            QNumbers.m(),
                            Theta,
                            Phi
                        );

                        PsiValue = Ylm * RadialPart;
                    }

                    Psi[{ ix, iz }] = PsiValue;
                }
            }

            Psi.normalize();
            return { Psi, QNumbers.hartreeEnergy() };
        }
    };
}