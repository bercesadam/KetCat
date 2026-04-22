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
#include "kwf_exporter/kwf_exporter.h"
#include "laser/laser_pulse.h"


using namespace KetCat;

int main() 
{
	constexpr natural_t NumBases = 5;
    constexpr natural_t NumQubits = 1;
    constexpr natural_t DiscretizationSteps = 512;
    constexpr real_t PhysicalExtent = 50.0;
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

    using Basis1D = basis_set_t<HilbertSpace1D, NumBases>;
    auto Bases = std::make_unique<Basis1D>(
    Basis1D  {{
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q0),
		EffectiveRadialOrbital<HilbertSpace1D, E>()(q1),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q2),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q3),
        EffectiveRadialOrbital<HilbertSpace1D, E>()(q4)
     }});

    auto DipoleMatrix = buildRadialDipoleMatrix(*Bases);

    for (natural_t i = 0; i < NumBases; ++i)
    {
        for (natural_t j = 0; j < NumBases; ++j)
        {
			std::cout << "Dipole Matrix Element μ(" << i << ", " << j << ") = " << DipoleMatrix[i][j].re << std::endl;
        }
	}


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
    auto Ortho = std::make_unique<Orthonormalizer<NumBases>>();
    Ortho->learn(*Bases2D);
    *Bases2D = Ortho->apply(*Bases2D);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    auto ReducedSpace =
        std::make_unique<ReducedEnergySpace<HilbertSpace, NumBases>>(*Bases2D);

    using ReducedHilbertSpace = typename decltype(ReducedSpace)::element_type::ReducedHilbertSpace;

    const natural_t Ket0Level = 0;
    const natural_t Ket1Level = 1;
    auto SeedReducedSpace = ReducedSpace->project((*Bases2D)[Ket0Level].m_Psi);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    const real_t dt = 0.01;
    real_t Time = 0.0;
    natural_t Frame = 0;

    const real_t HartreeEnergyDiff = (*Bases2D)[3].m_Energy - (*Bases2D)[2].m_Energy;
    std::cout << "Hartree energy difference between |0> and |1>: " << HartreeEnergyDiff << std::endl;
    //RwaRabiHamiltonian<NumBases> H(ReducedSpace->getEnergies(), DipoleMatrix, HartreeEnergyDiff, 0.01, 0.0);
	

    auto H = LaserHamiltonianBuilder<NumBases>::build(
        ReducedSpace->getEnergies(),
        DipoleMatrix,
        514.8,     // Wavelength in nm
        5e12,       // Intensity in W/cm²
        (*Bases2D)[0].m_Energy  // Reference level for rotating frame
	);
    tridiagonal_matrix_t<NumBases> Hmat = H.getMatrix();

    auto H2 = LaserHamiltonianBuilder<NumBases>::build(
        ReducedSpace->getEnergies(),
        DipoleMatrix,
        4500.0,     // Wavelength in nm
        5e12,       // Intensity in W/cm²
        (*Bases2D)[2].m_Energy  // Reference level for rotating frame
	);
    tridiagonal_matrix_t<NumBases> Hmat2 = H2.getMatrix();
    
    QuditSubspaceHelper<NumBases, NumQubits> Register;
    auto Psi = Register.productStateFromSeed(SeedReducedSpace);

    using OpSpace = decltype(Register)::OperationSpace<1>;
    using Solver = CrankNicolsonSolver<OpSpace>;
    Solver solver(Hmat, dt);

    StateVectorExporter<HilbertSpace> Exporter
    (
        "simulation.kwf",
        KetCat::ExportMode::RealImag
    );

    real_t LaserNm = 514.8;
    natural_t Skip = 10;
    while (Frame < 2000000)
    {
        //Psi = solver(Psi);

       Register.applyHamiltonian<1>(solver, Psi, {0}, Hmat);
       //Register.applyHamiltonian<1>(solver, Psi, {1}, Hmat);

        auto Psi0 = Register.extractLocalState(Psi, 0);
        //auto Psi1 = Register.extractLocalState(Psi, 1);
        
        std::ostringstream Title;
        Title << std::fixed << std::setprecision(2);
        Title << "Cesium (Z=55)|";
		Title << "Laser: " << LaserNm << " nm; 5E12 W/cm²|";
        Title << "Time: " << Time << " a.u.|";
        Title << "Populations: ";
		Title << "6s: " << Psi0[0].normSquared() * 100.0 << "% ";
		Title << "7p: " << Psi0[1].normSquared() * 100.0 << "% ";
		Title << "8d: " << Psi0[2].normSquared() * 100.0 << "% ";
		Title << "9f: " << Psi0[3].normSquared() * 100.0 << "% ";
		Title << "10g: " << Psi0[4].normSquared() * 100.0 << "%";
        /*
        if (Psi0[2].normSquared() >= 0.7)
        {
            Hmat = Hmat2;
            LaserNm = 4500.0;
            Skip = 1000;
            solver = Solver(Hmat, dt);
        }*/

        if (Frame % Skip == 0)
        {
            std::cout << Frame << std::endl;

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

            auto Psi_ = ReducedSpace->embed(Psi0);
            Exporter.writeTimestep(Time, Psi_, Title.str());
         }
        
        Time += dt;
        Frame++;
    }
}
