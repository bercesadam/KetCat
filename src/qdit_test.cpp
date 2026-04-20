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
#include "laser/laser_pulse.h"


using namespace KetCat;

int main() 
{
	constexpr natural_t NumBases = 5;
    constexpr natural_t NumQubits = 1;
    constexpr natural_t DiscretizationSteps = 256;
    constexpr real_t PhysicalExtent = 150.0;
    using HilbertSpace = InfiniteHilbertSpace<2_D, DiscretizationSteps, PhysicalExtent>;
    using HilbertSpace1D = InfiniteHilbertSpace<1_D, DiscretizationSteps, PhysicalExtent>;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    constexpr Element E = Element::Cs;

    using namespace SpectroscopicLetters;
    constexpr auto q0 = QuantumNumber<6, s>();
    constexpr auto q1 = QuantumNumber<7, p>();
    constexpr auto q2 = QuantumNumber<8, d>();
    constexpr auto q3 = QuantumNumber<9, f>();
    constexpr auto q4 = QuantumNumber<10, g>();

    basis_set_t<HilbertSpace1D, NumBases> Bases =
    { {
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q0),
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q1),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q2),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q3),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q4)
     } };

    auto DipoleMatrix = buildRadialDipoleMatrix(Bases);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    using Basis2D = basis_set_t<HilbertSpace, NumBases>;
    
    auto Bases2D = std::make_unique<Basis2D>(
        Basis2D{{
            Hydrogenic2D<HilbertSpace, E>()(q0),
            Hydrogenic2D<HilbertSpace, E>()(q1),
            Hydrogenic2D<HilbertSpace, E>()(q2),
            Hydrogenic2D<HilbertSpace, E>()(q3),
            Hydrogenic2D<HilbertSpace, E>()(q4)
        }}
    );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    auto ReducedSpace =
        std::make_unique<ReducedEnergySpace<HilbertSpace, NumBases>>(*Bases2D);

    using ReducedHilbertSpace = typename decltype(ReducedSpace)::element_type::ReducedHilbertSpace;

    const natural_t Ket0Level = 0;
    const natural_t Ket1Level = 1;
    auto SeedReducedSpace = ReducedSpace->project((*Bases2D)[Ket0Level].m_Psi);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    const real_t dt = 0.01;
    natural_t Frame = 0;

    const real_t HartreeEnergyDiff = (*Bases2D)[Ket1Level].m_Energy - (*Bases2D)[Ket0Level].m_Energy;
    std::cout << "Hartree energy difference between |0> and |1>: " << HartreeEnergyDiff << std::endl;
    //RwaRabiHamiltonian<NumBases> H(ReducedSpace->getEnergies(), DipoleMatrix, HartreeEnergyDiff, 0.01, 0.0);
	

    auto H = LaserHamiltonianBuilder<NumBases>::build(
        ReducedSpace->getEnergies(),
        DipoleMatrix,
        852.0,     // Wavelength in nm (Cs D2 line)
        1e6,       // Intensity in W/cm²
        Ket0Level  // Reference level for rotating frame
	);
    tridiagonal_matrix_t<NumBases> Hmat = H.getMatrix();
    
    QuditSubspaceHelper<NumBases, NumQubits> Register;
    auto Psi = Register.productStateFromSeed(SeedReducedSpace);

    using OpSpace = decltype(Register)::OperationSpace<1>;
    using Solver = CrankNicolsonSolver<OpSpace>;
    Solver solver(Hmat, dt);

    KetCat::StateVectorCsvExporter<HilbertSpace> Exporter
    (
        "simulation.csv",
        KetCat::ExportMode::RealImag
    );

    while (Frame < 35000)
    {
        //Psi = solver(Psi);

       Register.applyHamiltonian<1>(solver, Psi, {0}, Hmat);
       //Register.applyHamiltonian<1>(solver, Psi, {1}, Hmat);

        auto Psi0 = Register.extractLocalState(Psi, 0);
        //auto Psi1 = Register.extractLocalState(Psi, 1);

        auto Psi_ = ReducedSpace->embed(Psi0); 

        std::ostringstream Title;
        Title << "Cesium (Z=55) | Populations [";
		Title << std::fixed << std::setprecision(2);
		Title << "6s: " << Psi0[0].normSquared() * 100.0 << "% ";
		Title << "7p: " << Psi0[1].normSquared() * 100.0 << "% ";
		Title << "8d: " << Psi0[2].normSquared() * 100.0 << "% ";
		Title << "9f: " << Psi0[3].normSquared() * 100.0 << "% ";
		Title << "10g: " << Psi0[4].normSquared() * 100.0 << "%]";

        if (Frame % 20 == 0)
        {
            // Print all reduced space probabilities
            for (natural_t i = 0; i < NumBases; ++i)
            {
                std::cout << "Probability of basis state Qbit 0 " << i << ": " << Psi0[i].normSquared() * 100.0 << "%" << std::endl;

            }/*
            for (natural_t i = 0; i < NumBases; ++i)
            {
                std::cout << "Probability of basis state QBit 1 " << i << ": " << Psi1[i].normSquared() * 100.0 << "%" << std::endl;

            }*/
            std::cout << "------------------------" << std::endl;

            Exporter.writeTimestep(Frame, Psi_, Title.str());
         }

        Frame++;
    }
}
