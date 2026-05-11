#include <tuple>
#include <string>
#include <functional>

#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief Demonstration of hydrogen-like atomic orbitals in a 1D particle in a box system.

int main()
{
	constexpr real_t BoxLength = 1.0;
	constexpr natural_t Steps = 96;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;

	// Hydrogen orbital constructor with fixed parameters
	constexpr auto hydrogenCtor = std::bind(HydrogenOrbital<HilbertSpace>(), std::placeholders::_1, 0.05);

	// List of hydrogen orbitals to simulate: (StateVector, Name, l)
	using namespace SpectroscopicLetters;
	std::array<std::tuple<StateVector<HilbertSpace>, std::string, natural_t>, 6> hydrogenOrbitals =
	{
		 std::make_tuple(hydrogenCtor(QuantumNumber<1, s>()), "1s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber<2, p>()), "2p", 1),
		 std::make_tuple(hydrogenCtor(QuantumNumber<3, d>()), "3d", 2),
		 std::make_tuple(hydrogenCtor(QuantumNumber<4, f>()), "4f", 3),
		 std::make_tuple(hydrogenCtor(QuantumNumber<5, g>()), "5g", 4),
		 std::make_tuple(hydrogenCtor(QuantumNumber<6, h>()), "6h", 5)
	};

	constexpr real_t TimeStep = 1E-4;
	constexpr KetCat::real_t mass = 1.0;

	for (const auto& orbital : hydrogenOrbitals)
	{
		std::cout << "Hydrogen orbital: " << std::get<1>(orbital) << "\n";

		auto potential = SoftCoulombRadialPotential(
			1.0,		// Effective nuclear charge Z_eff
			2e-2,		// Softening parameter a
			std::get<2>(orbital), // Orbital quantum number ℓ
			hBar,		// Reduced Planck's constant ℏ
			mass	    // Reduced mass μ
		);

		auto hamiltonian = Hamiltonian<Steps>(mass, HilbertSpace::dx, potential);
		auto box = OneDimensionalParticleBox<HilbertSpace>(hamiltonian, std::get<0>(orbital), TimeStep);
		
		auto visu = Visu::VisuOscilloscope<Steps>(
			Visu::UsePhaseEncoding::NO,
			Visu::ClearScreen::NO,
			Visu::ShowComplexParts::NO,
			Visu::ShowPotential::YES
		);
		visu.setPotential(potential, HilbertSpace::dx);
		visu.update(box.evolve());
	}

	return 0;
}
