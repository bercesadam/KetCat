#pragma once
#include "hilbert_space/state_vector.h"
#include "hilbert_space/hilbert2d.h"

namespace KetCat
{

    template<dimension_t Dim, real_t Extent>
    struct Harmonic2D
    {
        using Space = Hilbert2D<Dim, Extent>;

        /**
         * @brief Generate 2D quantum harmonic oscillator seed wavefunction
         *
         * @param nx        Quantum number along x
         * @param ny        Quantum number along y
         *
         * @return StateVector<Space>
         */
        StateVector<Space> operator()(dimension_t nx, dimension_t ny)
        {
            StateVector<Space> Psi{ cplx_t::zero() };
            const real_t dx = Extent / (Dim - 1);

            auto Hermite = [](dimension_t n, real_t x) -> real_t
            {
                if (n == 0) return 1.0;
                if (n == 1) return 2.0 * x;
                real_t H0 = 1.0;
                real_t H1 = 2.0 * x;
                real_t Hn = 0.0;
                for (dimension_t i = 2; i <= n; ++i)
                {
                    Hn = 2.0 * x * H1 - 2.0 * (i - 1) * H0;
                    H0 = H1;
                    H1 = Hn;
                }
                return Hn;
            };

            // Simple alpha = 1 for normalized seed
            const real_t alpha = 1.0;

            for (dimension_t ix = 0; ix < Dim; ++ix)
            {
                for (dimension_t iy = 0; iy < Dim; ++iy)
                {
                    real_t x = (static_cast<real_t>(ix) - Dim / 2) * dx;
                    real_t y = (static_cast<real_t>(iy) - Dim / 2) * dx;

                    real_t Hx = Hermite(nx, alpha * x);
                    real_t Hy = Hermite(ny, alpha * y);
                    real_t envelope = ConstexprMath::exp<20>(-0.5 * alpha * alpha * (x * x + y * y));

                    cplx_t value = cplx_t(Hx * Hy * envelope, 0.0);

                    dimension_t index = Space::getIndex(ix, iy);
                    Psi[index] = value;
                }
            }

            return Psi;
        }
    };
}