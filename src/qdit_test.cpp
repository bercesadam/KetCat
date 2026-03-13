#include <SDL2/SDL.h>
#include <iostream>
#include <iomanip>
#include <memory>

#include "wavefunction/eigenstate.h"
//#include "wavefunction/2d_hydrogen.h"
#include "hilbert_space/reduced_energy_space.h"
#include "hamiltonian/qdit_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "visu/visu_oscilloscope.h"
#include "systems/qdit_register.h"

using namespace KetCat;
using namespace KetCat::Visu;

int main(int, char**) 
{
    constexpr natural_t N = 96; // Grid size
    constexpr real_t L = 10.0; // Physical size of the system
    constexpr natural_t NumBasis = 8; // Number of basis states for reduced space

    using HilbertSpace = InfiniteHilbertSpace<1_D, N, L>;
	EigenState<HilbertSpace> Generator;

    std::array<std::tuple<natural_t>, NumBasis> ParamSets;
    for (natural_t n = 0; n < NumBasis; ++n)
    {
       ParamSets[n] = std::make_tuple(n + 1);
    }

    /*
    Hydrogen2D<InfiniteHilbertSpace<2_D, N, L>> Generator;

    std::array<std::tuple<QuantumNumber>, NumBasis> ParamSets = {
        std::make_tuple(QuantumNumber::_1s()),
        std::make_tuple(QuantumNumber::_2p()),
        std::make_tuple(QuantumNumber::_3d()),
        std::make_tuple(QuantumNumber::_4f()),
        std::make_tuple(QuantumNumber::_5g()),
        std::make_tuple(QuantumNumber::_6h()),
        std::make_tuple(QuantumNumber::_7i()),
        std::make_tuple(QuantumNumber::_8k()),
    };
    */

    auto ReducedSpace =
        std::make_unique<ReducedEnergySpace<HilbertSpace, NumBasis>>(Generator, ParamSets);

    using ReducedHilbertSpace = typename decltype(ReducedSpace)::element_type::ReducedHilbertSpace;

    const natural_t Ket0Level = 0;
    const natural_t Ket1Level = 1;

    auto Wavefunction0 = Generator(std::get<0>(ParamSets[Ket0Level]));
    auto Wavefunction1 = Generator(std::get<0>(ParamSets[Ket1Level]));

    QuditSubspaceEngine<NumBasis, 3> Engine;
   


    // Initial state (|0>)
    auto Psi = ReducedSpace->project(Wavefunction0.m_Psi);

    using H = SingleQditGateHamiltonian<NumBasis>;
    const real_t HartreeEnergyDiff = Wavefunction1.m_Energy - Wavefunction0.m_Energy;
    std::cout << "Hartree energy difference between |0> and |1>: " << HartreeEnergyDiff << std::endl;
    const real_t Omega = 5.0;
    const real_t lambda = 0.1;
    const real_t leakageFactor = 0.0;
    auto Hx = H(RotationAxis::X, ReducedSpace->getEnergies(), Omega, lambda, leakageFactor, Ket0Level, Ket1Level);
    auto Hz = H(RotationAxis::Z, ReducedSpace->getEnergies(), Omega, lambda, leakageFactor, Ket0Level, Ket1Level);

    const real_t dt = 0.01;
    real_t time = 0.0;
    using Solver = CrankNicolsonSolver<ReducedHilbertSpace>;

    Visu::VisuOscilloscope<HilbertSpace::Dim> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::YES,
		Visu::ShowPotential::NO
	);

   // ===== Hadamard gate sequence =====

    // 1️⃣ X π/2 pulse
    const real_t tX = ConstexprMath::Pi / (Omega);
    real_t elapsed = 0.0;

    while (elapsed < tX)
    {
        tridiagonal_matrix_t Hmat = Hx(time);
        //Psi = Solver(Hmat, dt)(Psi);

        Engine.applyHamiltonian(Psi, {0}, Hmat, dt);

        auto Psi_ = ReducedSpace->embed(Psi);
        Psi_.normalize();

        std::cout << std::endl << "\nX |0> " << Psi_.probabilityOf(Wavefunction0.m_Psi)
                << "  |1> " << Psi_.probabilityOf(Wavefunction1.m_Psi)
                << "  Time: " << time << std::endl;

  //      visu.update(Psi_);

        time += dt;
        elapsed += dt;
    }
exit(0);
    // 2️⃣ Z π phase
    const real_t omega = HartreeEnergyDiff;
    const real_t tZ = ConstexprMath::Pi / omega;
    elapsed = 0.0;

    while (elapsed < tZ)
    {
        tridiagonal_matrix_t Hmat = Hz(time);
        Psi = Solver(Hmat, dt)(Psi);

        auto Psi_ = ReducedSpace->embed(Psi);
        Psi_.normalize();

        std::cout << std::endl << "\nZ |0> " << Psi_.probabilityOf(Wavefunction0.m_Psi)
                << "  |1> " << Psi_.probabilityOf(Wavefunction1.m_Psi)
                << "  Time: " << time << std::endl;

        //visu.update(Psi_);

        time += dt;
        elapsed += dt;
    }
}

