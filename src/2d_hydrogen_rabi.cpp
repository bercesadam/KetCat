#include <SDL2/SDL.h>
#include <iostream>
#include <iomanip>

#include "wavefunction/2d_hydrogen.h"
#include "visu/wavefunction_viewer/wavefunction_viewer.h"

using namespace KetCat;
using namespace KetCat::Visu;

/// This demo program is a smoke test of new 2D functionality.
/// It simulates and visualizes Rabi oscillations between two hydrogenic eigenstates (|3p₋₁⟩ and |4d₀⟩) in a 2D Hilbert space. 

int main(int, char**) 
{
    // Construct Hilbert-space
    constexpr natural_t DiscretizationSteps = 256;
    constexpr real_t PhysicalExtent = 20.0;
    using HilbertSpace = InfiniteHilbertSpace<2_D, DiscretizationSteps, PhysicalExtent>;

    // Select quantum numbers:
    // |3p, m = -1⟩ and |4d, m = 0⟩
    constexpr auto q0 = QuantumNumber::_3p_m1();
    constexpr auto q1 = QuantumNumber::_4d_m0();

    // Construct initial eigenstates:
    auto Psi0 = Hydrogen2D<HilbertSpace>()(q0);
    auto Psi1 = Hydrogen2D<HilbertSpace>()(q1);

    // Corresponding energy eigenvalues (Hartree units):
    const real_t E0 = q0.hartreeEnergy();
    const real_t E1 = q1.hartreeEnergy();

    // Rabi frequency Ω
    const real_t Omega = 0.02;

    // Time step Δt
    const real_t dt = 0.5;

    // Initialize visualization window
    bool Running = true;
    SDL_Event Event;
    WavefunctionViewer<HilbertSpace> Visu(1200, 900);
    int Frame = 0;

    // Initialize probabilities and the state vector
    real_t Palpha = 0.0; // Occupation probability P(t) = |α|²
    real_t Pbeta = 0.0;  // Occupation probability P(t) = |β|²
    auto Psi = Psi0;     // Initial state: pure ψ₀

    while (Running)
    {
        Frame++;
        real_t t = Frame * dt;   // Physical time

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
        cplx_t Phase0 = ConstexprMath::exp<20>(cplx_t(0.0, -E0 * t));
        cplx_t Phase1 = ConstexprMath::exp<20>(cplx_t(0.0, -E1 * t));

        // Rabi mixing angle:
        //
        // θ(t) = Ω t
        //
        // Governs population transfer between states
        double Theta = Omega * t;

        // Superposition coefficients:
        //
        // α(t) = e^{-i E₀ t} cos(θ/2)
        // β(t) = -i e^{-i E₁ t} sin(θ/2)
        //
        // The factor -i ensures proper phase relation
        // for unitary two-level rotation.
        cplx_t Alpha = Phase0 * cplx_t(std::cos(Theta / 2.0), 0.0);
        cplx_t Beta = Phase1 * cplx_t(0.0, -std::sin(Theta / 2.0));

        // While transition is not complete:
        if (Pbeta < 0.99)
        {
            // Construct superposition:
            //
            // ψ = α ψ₀ + β ψ₁
            Psi = Psi0.superpose(Psi1, Alpha, Beta);

            // Renormalize in discrete 2D space:
            //
            // ∑ |ψ|² dx² = 1
            Psi.normalize();
        }

        // Calculate occupation probabilities of state ψ₀ and ψ₁:
        Palpha = Psi.probabilityOf(Psi0);
        Pbeta = Psi.probabilityOf(Psi1);

        // Create custom title for the Visu
        std::ostringstream ProbaPopulation;
        ProbaPopulation << std::setprecision(2)
            << "Population: |⟨4d₀|ψ⟩|² = "
            << Pbeta * 100.0 << "%, "
            << "|⟨3p₁|ψ⟩|² = "
            << Palpha * 100.0 << "%";

        // Render density |ψ(x,y,t)|² and state coefficients
        Visu.render(Psi, ProbaPopulation.str());

        while (SDL_PollEvent(&Event))
        {
            if (Event.type == SDL_QUIT)
            {
                Running = false;
            }
        }
    }

    return 0;
}
