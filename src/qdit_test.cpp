#include <iostream>
#include <iomanip>
#include <memory>

#include "wavefunction/atomic/2d_hydrogenic.h"
#include "hilbert_space/reduced_energy_space.h"
#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "visu/visu_oscilloscope.h"
#include "systems/qdit_register.h"
#include "hilbert_space/gram_schmidt_orthonorm.h"
#include "hilbert_space/dipole_operator.h"
#include "visu/file_exporter.h"


using namespace KetCat;

int main() 
{
	constexpr natural_t NumBases = 5;
    constexpr natural_t DiscretizationSteps = 256;
    constexpr real_t PhysicalExtent = 150.0;
    using HilbertSpace = InfiniteHilbertSpace<2_D, DiscretizationSteps, PhysicalExtent>;
    using HilbertSpace1D = InfiniteHilbertSpace<1_D, DiscretizationSteps, PhysicalExtent>;

    constexpr Element E = Element::H;

    using namespace SpectroscopicLetters;
    constexpr auto q0 = QuantumNumber<1, s, 0>();
    constexpr auto q1 = QuantumNumber<2, p, 0>();
    constexpr auto q2 = QuantumNumber<3, d, 0>();
    constexpr auto q3 = QuantumNumber<4, f, 0>();
    constexpr auto q4 = QuantumNumber<5, g, 0>();

    constexpr basis_set_t<HilbertSpace1D, NumBases> Bases =
    { {
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q0),
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q1),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q2),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q3),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q4)
     } };

    auto DipoleMatrix = buildRadialDipoleMatrix(Bases);

    basis_set_t< HilbertSpace, NumBases> Bases2D =
    { {
        Hydrogenic2D<HilbertSpace, E>()(q0),
        Hydrogenic2D<HilbertSpace, E>()(q1),
        Hydrogenic2D<HilbertSpace, E>()(q2),
        Hydrogenic2D<HilbertSpace, E>()(q3),
        Hydrogenic2D<HilbertSpace, E>()(q4)
     } };

	std::cout << "Dipole matrix elements (x direction):" << std::endl;
    for (natural_t n = 0; n < Bases.size(); ++n)
    {
        for (natural_t m = 0; m < Bases.size(); ++m)
        {
            std::cout << std::setw(20) << n << ";" << m << ":" << DipoleMatrix[n][m].re << " ";
        }
        std::cout << std::endl;
	}

    auto ReducedSpace =
        std::make_unique<ReducedEnergySpace<HilbertSpace, NumBases>>(Bases2D);

    using ReducedHilbertSpace = typename decltype(ReducedSpace)::element_type::ReducedHilbertSpace;

    const natural_t Ket0Level = 0;
    const natural_t Ket1Level = 1;

    auto const Wavefunction0 = Bases2D[Ket0Level];
    auto const Wavefunction1 = Bases2D[Ket1Level];

    //QuditSubspaceHelper<NumBasis, 3> Helper;
    
    // Initial state (|0>)
    auto Seed = ReducedSpace->project(Wavefunction0.m_Psi);
    auto Psi = Seed; 

    const real_t HartreeEnergyDiff = Wavefunction1.m_Energy - Wavefunction0.m_Energy;
    std::cout << "Hartree energy difference between |0> and |1>: " << HartreeEnergyDiff << std::endl;

    const real_t dt = 0.01;
    natural_t Frame = 0;

    RwaRabiHamiltonian<NumBases> H(ReducedSpace->getEnergies(), DipoleMatrix, HartreeEnergyDiff, 0.01, 0.0);
	tridiagonal_matrix_t<NumBases> Hmat = H.getMatrix();
    using Solver = CrankNicolsonSolver<ReducedHilbertSpace>;
    Solver solver(Hmat, dt);

    /*
    Visu::VisuOscilloscope<HilbertSpace::Dim> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::YES,
		Visu::ShowPotential::NO
	);*/

    KetCat::StateVectorCsvExporter<HilbertSpace> exporter
    (
        "simulation.csv",
        KetCat::ExportMode::RealImag
    );

    while (true)
    {
        Psi = solver(Psi);

        //Helper.applyHamiltonian<1>(Psi, {0}, Hmat, dt);

        auto Psi_ = ReducedSpace->embed(Psi); //Helper.extractLocalState(Psi, 0));

        //Psi_.normalize();

        if (Frame % 100 == 0)
        {
            // Print all reduced space probabilities
            for (natural_t i = 0; i < NumBases; ++i)
            {
                std::cout << "Probability of basis state " << i << ": " << Psi[i].normSquared() * 100.0 << "%" << std::endl;

            }
            std::cout << "------------------------" << std::endl;
         }

        if (Frame % 100 == 0) exporter.writeTimestep(Frame, Psi_, "test");
        Frame++;
    }
}
