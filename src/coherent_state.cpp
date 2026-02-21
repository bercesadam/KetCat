#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@file Demonstration of a coherent state in a 1D harmonic oscillator potential, probably the most
///      classical-like quantum state possible and in the same time the most impressive visualization
///      effects which can be witnessed in this simple simulator project.
///      Optionally, the width σ of the Gaussian can be overridden to observe "quantum breathing" effects.

int main()
{
	constexpr real_t BoxLength = 3.0;
	constexpr natural_t Steps = 96;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;


	constexpr KetCat::real_t center = BoxLength / 2.0;

	constexpr auto psi0 =
		CoherentStateGaussian<HilbertSpace>()(
			KetCat::hBar,	// Reduced Planck constant
			center,			// x₀
			10.0,			// p₀
			1.0,			// m
			1.0				// ω
//			1.4 			// σ override - optional, for "quantum breathing" effect
		);

	constexpr auto potential = HarmonicOscillatorPotential(
		1.0,    // m
		1.0,    // ω
		center  // x₀
	);

	constexpr real_t TimeStep = 5E-4;
	constexpr KetCat::real_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<HilbertSpace::Dim>(mass, HilbertSpace::dx, potential);

	OneDimensionalParticleBox<HilbertSpace> box(hamiltonian, psi0, TimeStep);
	
	Visu::VisuOscilloscope<Steps> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::YES,
		Visu::ShowPotential::YES
	);

	visu.setPotential(potential, HilbertSpace::dx);

	while (true)
	{
		auto p = box.evolve();
		p.normalize();
		visu.update(p);
	}

	return 0;
}
