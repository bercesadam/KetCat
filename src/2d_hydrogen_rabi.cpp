#include <iostream>
#include <iomanip>
#include <memory>

#include "wavefunction/atomic/2d_hydrogenic.h"
#include "visu/file_exporter.h"

using namespace KetCat;

/// This demo program is a smoke test of new 2D functionality.
/// It simulates and visualizes Rabi oscillations between two hydrogenic eigenstates (|3p₋₁⟩ and |4d₀⟩) in a 2D Hilbert space. 

int main() 
{
    // Construct Hilbert-space
    constexpr natural_t DiscretizationSteps = 256;
    constexpr real_t PhysicalExtent = 200.0;
    using HilbertSpace = InfiniteHilbertSpace<2_D, DiscretizationSteps, PhysicalExtent, GridType::Logarithmic>;

    // Select quantum numbers:
    // |3p, m = -1⟩ and |4d, m = 0⟩
    using namespace SpectroscopicLetters;
	constexpr auto q0 = QuantumNumber<6, s, 0>();
	constexpr auto q1 = QuantumNumber<7, p, 0>();

    // Construct initial eigenstates:
    auto Psi0 = std::make_unique<StateVector<HilbertSpace>>(
        Hydrogen2D<HilbertSpace, Element::Cs>()(q0).m_Psi
    );

    auto Psi1 = std::make_unique<StateVector<HilbertSpace>>(
        Hydrogen2D<HilbertSpace, Element::Cs>()(q1).m_Psi
    );


    std::cout << "Psi0 self-overlap: " << Psi0->innerProduct(*Psi0).re << "\n";
std::cout << "Psi1 self-overlap: " << Psi1->innerProduct(*Psi1).re << "\n";
std::cout << "<Psi0|Psi1>: "       << Psi0->innerProduct(*Psi1).re << "\n";

    // Corresponding energy eigenvalues (Hartree units):
    const real_t E0 = q0.hartreeEnergy();
    const real_t E1 = q1.hartreeEnergy();

    // Rabi frequency Ω
    const real_t Omega = 0.02;

    // Time step Δt
    const real_t dt = 0.5;

    int Frame = 0;

    // Initialize probabilities and the state vector
    real_t Palpha = 0.0; // Occupation probability P(t) = |α|²
    real_t Pbeta = 0.0;  // Occupation probability P(t) = |β|²

    KetCat::StateVectorCsvExporter<HilbertSpace> exporter(
        "simulation.csv",
        KetCat::ExportMode::RealImag
    );

    while (Frame < 300)
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
        complex_t Phase0 = ConstexprMath::exp(complex_t(0.0, -E0 * t));
        complex_t Phase1 = ConstexprMath::exp(complex_t(0.0, -E1 * t));

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
        complex_t Alpha = Phase0 * complex_t(std::cos(Theta / 2.0), 0.0);
        complex_t Beta = Phase1 * complex_t(0.0, -std::sin(Theta / 2.0));

        // While transition is not complete:
        if (Pbeta < 0.99)
        {
            // Construct superposition:
            //
            // ψ = α ψ₀ + β ψ₁
            auto Psi = Psi0->superpose(*Psi1, Alpha, Beta);

            // Renormalize in discrete 2D space:
            //
            // ∑ |ψ|² dx² = 1
            Psi.normalize();

            // Calculate occupation probabilities of state ψ₀ and ψ₁:
            Palpha = Psi.probabilityOf(*Psi0);
            Pbeta = Psi.probabilityOf(*Psi1);

            // Create custom title for the Visu
            std::ostringstream ProbaPopulation;
            ProbaPopulation << std::setprecision(2)
                << "Population: |⟨4d₀|ψ⟩|² = "
                << Pbeta * 100.0 << "%, "
                << "|⟨3p₁|ψ⟩|² = "
                << Palpha * 100.0 << "%";
			std::cout << ProbaPopulation.str() << std::endl;

            // Render density |ψ(x,y,t)|² and state coefficients
            exporter.writeTimestep(t, Psi);
        }
    }

    return 0;
}
