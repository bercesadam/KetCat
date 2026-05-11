#pragma once
#include "wavefunction/wavefunction.h"
#include "effective_radial_orbital.h"

#include "atomic_physics_core/quantum_number.h"
#include "atomic_physics_core/elements.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

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
                ConstexprMath::factorial(L - ConstexprMath::abs(M)) /
                ConstexprMath::factorial(L + ConstexprMath::abs(M)) *
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
    inline complex_t sphericalHarmonic(int L, int M, real_t Theta, real_t Phi)
    {
        real_t Norm =
            ConstexprMath::sqrt((2.0 * L + 1.0) /
                      (4.0 * ConstexprMath::Pi) *
                      ConstexprMath::factorial(L - ConstexprMath::abs(M)) /
                      ConstexprMath::factorial(L + ConstexprMath::abs(M)));

        real_t X = ConstexprMath::cos(Theta);
        real_t P = legendre(L, ConstexprMath::abs(M), X);

        // Phase factor e^{i m φ}
        complex_t Phase = ConstexprMath::exp(complex_t(0.0, M * Phi));

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
    template<spatial_hilbert_space_with_dim_t<2_D> HilbertSpace, Element element>
    struct Hydrogenic2D
    {
        
        /// @brief Generate a 2D hydrogenic orbital for (n,l,m).
        /// @param QNumbers QuantumNumber object containing (n,l,m).
        /// @return StateVector<HilbertSpace> containing Ψ(x,z) values.
		template <quantum_number_t QuantumNumberType>
        constexpr Wavefunction<HilbertSpace> operator()(QuantumNumberType QNumbers) noexcept
        {
            constexpr natural_t Steps = HilbertSpace::Steps;

            StateVector<HilbertSpace> Psi{ complex_t::zero() };

            constexpr real_t QuantumDefect = RydbergQuantumDefect::value(element, QNumbers);
            constexpr real_t EffectiveN = static_cast<real_t>(QNumbers.n()) - QuantumDefect;

            // Experimental feature: scale the radial part based on the distance from the nucleus,
			// so Rydberg states can fit in the same grid as low-lying states, while still showing the correct nodal structure.
            constexpr real_t Scale = EffectiveN * EffectiveN;
            constexpr real_t A = 80.0;
			constexpr real_t B = 0.2;
			constexpr real_t C = 2;
            constexpr real_t Factor = ((A / (Scale + B)) + C);
            constexpr real_t PhysicalExtent = Scale * Factor;

            // --------------------------------------------------------
            // 1D radial component Rₙₗ(r) = uₙₗ(r) / r
            // Already normalized by HydrogenOrbital<Dim>().
            // --------------------------------------------------------
            using RadialSpace = InfiniteHilbertSpace<1_D, Steps, HilbertSpace::Extent, HilbertSpace::Grid>;
            auto RadialWavefunction = EffectiveRadialOrbital<RadialSpace, element>{}(QNumbers);
            auto RadialArray = RadialWavefunction.m_Psi;

            const real_t ViewScale = 1.0; //PhysicalExtent / HilbertSpace::Extent;

			std::cout << Factor << ", " << PhysicalExtent << std::endl;
            for (natural_t ix = 0; ix < Steps; ++ix)
            {
                for (natural_t iz = 0; iz < Steps; ++iz)
                {
                    // --------------------------------------------
                    // Physical grid coordinates (centered)
                    // --------------------------------------------
                    natural_t ixFromCenter = (ix < Steps / 2)
                        ? (Steps / 2 - ix)
                        : (ix - Steps / 2);

                    natural_t izFromCenter = (iz < Steps / 2)
                        ? (Steps / 2 - iz)
                        : (iz - Steps / 2);

                    real_t XCoord = ((ix < Steps / 2) ?
                        -HilbertSpace::gridToR(ixFromCenter) :
                         HilbertSpace::gridToR(ixFromCenter))
                        * ViewScale;

                    real_t ZCoord = ((iz < Steps / 2) ?
                        -HilbertSpace::gridToR(izFromCenter) :
                         HilbertSpace::gridToR(izFromCenter))
                        * ViewScale;

                    real_t YCoord = HilbertSpace::RMin * ViewScale;

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
                        natural_t IR = RadialSpace::rToGrid(R);
                        
                        // Zero-clamp the radial overflow 
                        if (IR >= Steps)
                        {
                            Psi[{ ix, iz }] = complex_t::zero();
                            continue;
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
                        complex_t Ylm = sphericalHarmonic(
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
            return { Psi, RadialWavefunction.m_Energy };
        }
    };
}