#include <tuple>
#include <string>
#include <functional>

#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief Demonstration of hydrogen-like atomic orbitals in a 1D particle in a box system.

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

	// Hydrogen orbital constructor with fixed parameters
	constexpr auto hydrogenCtor = std::bind(HydrogenOrbital<cfg.M>(), std::placeholders::_1, 0.05, cfg.dx);

	// List of hydrogen orbitals to simulate: (StateVector, Name, l)
	std::array<std::tuple<StateVector<InfiniteHilbertSpace<94>>, std::string, unsigned int>, 6> hydrogenOrbitals =
	{
		 std::make_tuple(hydrogenCtor(QuantumNumber::_1s()), "1s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2s()), "2s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2p()), "2p", 1),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3s()), "3s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3p()), "3p", 1),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3d()), "3d", 2)
	};

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

		auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potential);

		auto box = OneDimensionalParticleBox<cfg.N>(
			cfg,
			hamiltonian,
			std::get<0>(orbital)
		);
		
		auto visu = Visu::VisuOscilloscope<cfg.M>(
			Visu::UsePhaseEncoding::NO,
			Visu::ClearScreen::NO,
			Visu::ShowComplexParts::NO,
			Visu::ShowPotential::YES
		);
		visu.setPotential(potential, cfg.dx);
		visu.update(box.evolve());
	}

	return 0;
}
