#include <SDL2/SDL.h>
#include <vector>
#include <iostream>
#include "visu/wavefunction_viewer/wavefunction_viewer.h"
#include "wavefunction/2d_hydrogen.h"
#include "wavefunction/2d_harmonic.h"

using namespace KetCat;
using namespace KetCat::Visu;

int main() 
{
    constexpr int N = 200;
    constexpr real_t ex = 20.0;
    using Space = Hilbert2D<N, ex>;

    // Select quantum numbers:
    // |3p, m = -1⟩ and |4d, m = 0⟩
    auto q0 = QuantumNumber::_3p_m1();
    auto q1 = QuantumNumber::_4d_m0();

    // Construct stationary eigenstates:
    // ψ₀(x,y), ψ₁(x,y)
    auto psi0 = Hydrogen2D<N, ex>()(q0);
    auto psi1 = Hydrogen2D<N, ex>()(q1);

    // Corresponding energy eigenvalues (Hartree units):
    // Eₙ from:  H ψₙ = Eₙ ψₙ
    const double E0 = q0.hartreeEnergy();
    const double E1 = q1.hartreeEnergy();

    // Rabi frequency Ω
    const double Omega = 0.02;

    // Time step Δt
    const double dt = 0.5;

    // Visualization window
    WavefunctionViewer<Space> visu(1200, 900);

    int frame = 0;
    double P = 0.0;     // Occupation probability P(t) = |β|²
    auto psi = psi0;    // Initial state: pure ψ₀

    while (true)
    {
        frame++;
        double t = frame * dt;   // Physical time

        // Time evolution phase factors:
        //
        // phase₀ = e^{-i E₀ t}
        // phase₁ = e^{-i E₁ t}
        //
        // Derived from Schrödinger equation:
        //   i ∂ψ/∂t = H ψ
        //
        // For energy eigenstates:
        //   ψ(t) = ψ(0) e^{-i E t}
        cplx_t phase0 = ConstexprMath::exp<20>(cplx_t(0.0, -E0 * t));
        cplx_t phase1 = ConstexprMath::exp<20>(cplx_t(0.0, -E1 * t));


        // Rabi mixing angle:
        //
        // θ(t) = Ω t
        //
        // Governs population transfer between states
        double theta = Omega * t;

        // Superposition coefficients:
        //
        // α(t) = e^{-i E₀ t} cos(θ/2)
        // β(t) = -i e^{-i E₁ t} sin(θ/2)
        //
        // The factor -i ensures proper phase relation
        // for unitary two-level rotation.
        cplx_t alpha = phase0 * cplx_t(std::cos(theta / 2.0), 0.0);
        cplx_t beta = phase1 * cplx_t(0.0, -std::sin(theta / 2.0));

        // While transition is not complete:
        if (P < 0.99)
        {
            // Construct superposition:
            //
            // ψ = α ψ₀ + β ψ₁
            psi = psi0.superpose(psi1, alpha, beta);

            // Renormalize in discrete 2D space:
            //
            // ∑ |ψ|² dx² = 1
            psi.normalize2D(Space::dx);
        }

        // Occupation probability of state ψ₁:
        //
        // P(t) = |β(t)|²
        //
        // Expected analytic behavior:
        //   P(t) = sin²(Ω t / 2)
        P = beta.normSquared();

        // Render density |ψ(x,y,t)|² and state coefficients
        visu.render(psi, alpha, beta);
    }

    return 0;
}
