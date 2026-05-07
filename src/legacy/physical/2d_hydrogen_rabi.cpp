#include <iostream>
#include <iomanip>
#include <memory>

#include "hilbert_space/gram_schmidt_orthonorm.h"
#include "wavefunction/atomic/2d_hydrogenic.h"
#include "visu/file_exporter.h"

using namespace KetCat;

/// Smoke test of hydrogenic wavefunctions (STO and QDT) in 2D, simulating faux Rabi oscillations between adjacent states

int main() 
{
    // Construct Hilbert-space
    constexpr natural_t DiscretizationSteps = 256;
    constexpr real_t PhysicalExtent = 200.0;
    using HilbertSpace = InfiniteHilbertSpace<2_D, DiscretizationSteps, PhysicalExtent>;

    using namespace SpectroscopicLetters;
	constexpr auto q0 = QuantumNumber<6, s, 0>();
	constexpr auto q1 = QuantumNumber<7, p, 0>();
    constexpr auto q2 = QuantumNumber<8, d, 0>();
    constexpr auto q4 = QuantumNumber<9, f, 0>();
    constexpr auto q5 = QuantumNumber<10, g, 0>();

    basis_set_t<HilbertSpace, 5> Seed =   
    {{
        Hydrogenic2D<HilbertSpace, Element::Cs>()(q0),
        Hydrogenic2D<HilbertSpace, Element::Cs>()(q1),
        Hydrogenic2D<HilbertSpace, Element::Cs>()(q2),
        Hydrogenic2D<HilbertSpace, Element::Cs>()(q4),
        Hydrogenic2D<HilbertSpace, Element::Cs>()(q5)
     }};

	Orthonormalizer<5, true> Ortho;
	Ortho.learn(Seed);
	auto Basis = Ortho.apply(Seed);

    std::array<std::string, 4> Captions =
    {
		"Cesium (Z=55) 6s (STO) -> 7p (STO)",
		"Cesium (Z=55) 7p (STO) -> 8d (STO)",
		"Cesium (Z=55) 8d (STO) -> 9f (QDT)",
		"Cesium (Z=55) 9f (STO) -> 10g (QDT)",
    };

    KetCat::StateVectorCsvExporter<HilbertSpace> exporter
    (
        "simulation.csv",
        KetCat::ExportMode::RealImag
    );

    // Rabi frequency Ω
    const real_t Omega = 0.02;

    // Time step Δt
    const real_t dt = 0.5;

    auto Psi0 = Basis[0].m_Psi;

    for (natural_t i = 0; i < Basis.size() - 1; ++i)
    {
        int Frame = 0;

        real_t E0 = Basis[i].m_Energy;
        real_t E1 = Basis[i + 1].m_Energy;
        
        auto Psi1 = Basis[i + 1].m_Psi;

        // Initialize probabilities and the state vector
        real_t Palpha = 0.0; // Occupation probability P(t) = |α|²
        real_t Pbeta = 0.0;  // Occupation probability P(t) = |β|²

        while (Pbeta < 0.999)
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

            // Construct superposition:
            //
            // ψ = α ψ₀ + β ψ₁
            auto Psi = Psi0.superpose(Psi1, Alpha, Beta);

            // Renormalize in discrete 2D space:
            //
            // ∑ |ψ|² dx² = 1
            Psi.normalize();

            // Calculate occupation probabilities of state ψ₀ and ψ₁:
            Palpha = Psi.probabilityOf(Psi0);
            Pbeta = Psi.probabilityOf(Psi1);

            std::cout << Captions[i] << std::endl;
            if (Frame % 3 == 0) exporter.writeTimestep(t, Psi, Captions[i]);
        }

		Psi0 = Psi1; // Prepare for next transition
    }

    return 0;
}
