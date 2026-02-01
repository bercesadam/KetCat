#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@file Demonstration of a coherent state in a 1D harmonic oscillator potential, probably the most
///      classical-like quantum state possible and in the same time the most impressive visualization
///      effects which can be witnessed in this simple simulator project.
///      Optionally, the width σ of the Gaussian can be overridden to observe "quantum breathing" effects.

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(3.0, 5E-4);

	constexpr KetCat::float_t center = cfg.L / 2.0;

	constexpr auto psi0 =
		CoherentStateGaussian<cfg.M>()(
			KetCat::hBar,	// Reduced Planck constant
			center,			// x₀
			10.0,			// p₀
			cfg.dx,			// dx
			1.0,			// m
			1.0				// ω
//			1.4 			// σ override - optional, for "quantum breathing" effect
		);

	constexpr auto potential = HarmonicOscillatorPotential(
		1.0,    // m
		1.0,    // ω
		center  // x₀
	);

	constexpr KetCat::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potential);

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		psi0
	);
	
	Visu::VisuOscilloscope<cfg.M> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::YES,
		Visu::ShowPotential::YES
	);

	while (true)
	{
		auto p = box.evolve();
		p.normalize_with_dx(cfg.dx);
		visu.update(p, potential, cfg.dx);
	}

	return 0;
}
