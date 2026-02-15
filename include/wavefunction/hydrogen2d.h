#pragma once
#include "hilbert_space/state_vector.h"
#include "hilbert_space/hilbert2d.h"
#include "quantum_number.h"
#include "hydrogen.h"

namespace KetCat
{
    real_t legendre(int l, int m, real_t x)
    {
        if (m < 0)
        {
            real_t sign = (m % 2 == 0) ? 1.0 : -1.0;
            return sign * ConstexprMath::factorial(l - std::abs(m)) /
                          ConstexprMath::factorial(l + std::abs(m)) * legendre(l, -m, x);
        }

        real_t p_mm = 1.0;
        if (m > 0) {
            real_t somx2 = ConstexprMath::sqrt((1.0 - x) * (1.0 + x));
            real_t fact = 1.0;
            for (int i = 1; i <= m; ++i) { p_mm *= -fact * somx2; fact += 2; }
        }

        if (l == m)  return p_mm;

        real_t p_m1m = x * (2 * m + 1) * p_mm;
        if (l == m + 1) return p_m1m;

        real_t p_lm_prev = p_mm;
        real_t p_lm = p_m1m;
        for (int ll = m + 2; ll <= l; ++ll) {
            real_t p_next = ((2 * ll - 1) * x * p_lm - (ll + m - 1) * p_lm_prev) / (ll - m);
            p_lm_prev = p_lm;
            p_lm = p_next;
        }
        return p_lm;
    };

    /// Gömbi harmonikus Y_l^m(θ,φ)
    cplx_t spherical_harmonic(int l, int m, real_t theta, real_t phi)
    {
        real_t norm = std::sqrt((2.0 * l + 1.0) / (4.0 * ConstexprMath::Pi)
            * ConstexprMath::factorial(l - ConstexprMath::abs(m)) / ConstexprMath::factorial(l + ConstexprMath::abs(m)));
        real_t x = ConstexprMath::cos(theta);
        real_t P = legendre(l, ConstexprMath::abs(m), x);
        cplx_t phase = ConstexprMath::exp<20>(cplx_t(0.0, m * phi));
        if (m < 0) return phase.conj() * (((m % 2 == 0) ? 1.0 : -1.0) * norm * P);
        return phase * (norm * P);
    }


    /**
        * ---------------------------------------------------------------------------
        * @brief Hydrogen atom 2D slice wavefunction generator
        * ---------------------------------------------------------------------------
        *
        * Generates a planar slice of the full 3D hydrogen wavefunction:
        *
        *      Ψₙₗₘ(r,θ,φ) = Rₙₗ(r) · Yₗᵐ(θ,φ)
        *
        * where:
        *
        *      r  = √(x² + y² + z₀²)
        *      θ  = arccos(z₀ / r)
        *      φ  = atan2(y, x)
        *
        * and:
        *
        *      Rₙₗ(r) = uₙₗ(r) / r
        *
        * The radial component uₙₗ(r) is obtained from:
        *
        *      HydrogenOrbital<Dim>()(q, a_eff, dx)
        *
        * ---------------------------------------------------------------------------
        * Geometry
        * ---------------------------------------------------------------------------
        *
        * The generated state corresponds to the fixed plane:
        *
        *      z = zSlice
        *
        * Typically:
        *
        *      zSlice = 0   → XY plane
        *
        * Produces nice orbital shapes:
        *
        *      p-orbitals  →  butterfly
        *      d-orbitals  →  cloverleaf
        *      f-orbitals  →  flower patterns
        */

    template<dimension_t Dim, real_t Extent>
    struct Hydrogen2D
    {
        using Space = Hilbert2D<Dim, Extent>;

        /**
            * @brief Generate 2D hydrogen orbital slice
            *
            * @param q         Quantum numbers (n,l,m)
            * @param zSlice    Fixed z-plane value
            *
            * @return StateVector<Space>
            */
        StateVector<Space> operator()(QuantumNumber q)
        {
            StateVector<Space> Psi{ cplx_t::zero() };

            const real_t dx = Extent / (Dim - 1);

            // --------------------------------------------------------
            // Radial 1D orbital (already normalized in your API)
            // --------------------------------------------------------
            auto radial = HydrogenOrbital<Dim>()(q, 1.0, dx);

            for (dimension_t ix = 0; ix < Dim; ++ix)
            {
                for (dimension_t iz = 0; iz < Dim; ++iz)
                {
                    real_t x = (static_cast<real_t>(ix) - Dim / 2) * dx;
                    real_t z = (static_cast<real_t>(iz) - Dim / 2) * dx;
                    real_t y = 0.01;

                    real_t r = ConstexprMath::sqrt(x * x + y * y + z * z);

                    cplx_t value = cplx_t::zero();

                    if (r > 1e-12)
                    {
                        // --------------------------------------------
                        // Radial lookup
                        // --------------------------------------------
                        dimension_t ir = static_cast<dimension_t>(r / dx);
                        if (ir >= Dim)
                            ir = Dim - 1;

                        real_t u = radial[ir].re;   // magnitude only
                        real_t R = u / r;

                        // --------------------------------------------
                        // Angular coordinates
                        // --------------------------------------------
                        real_t theta = ConstexprMath::acos(z / r);
                        real_t phi = ConstexprMath::atan2(y, x);

                        cplx_t Y = spherical_harmonic(q.l(), q.m(), theta, phi);

                        value = Y * R;
                    }

                    dimension_t index = Space::getIndex(ix, iz);
                    Psi[index] = value;
                }
            }

            return Psi;
        }
    };

}
