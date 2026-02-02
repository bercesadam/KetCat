#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief A classical quantum physics textbook example: a Gaussian wave packet in an infinite potential well
/// with a barrier in the middle (can be configured or removed to see free propagation).
/// This example demonstrates quantum tunneling through the potential barrier, which is generally the way
/// FGMOS components in SSDs work in practice (electrons tunneling through oxide barriers).

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

    constexpr KetCat::float_t x0 = 0.01;
	constexpr KetCat::float_t sigma = 0.1;
	constexpr KetCat::float_t k0 = ConstexprMath::Pi * 10;

	constexpr auto gaussianPacKetCat = FreeParticleGaussianWavePacket<cfg.M>()(x0, k0, sigma, cfg.dx);

	constexpr PotentialBarrier potentialBarrier{
		0.45, 0.55, // potential wall in the middle
		3000
	};

	constexpr KetCat::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potentialBarrier);

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		gaussianPacKetCat
	);
	
	Visu::VisuOscilloscope<cfg.M> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::NO,
		Visu::ShowPotential::YES
	);

	visu.setPotential(potentialBarrier, cfg.dx);

	while (true)
	{
		auto p = box.evolve();
		p.normalize_with_dx(cfg.dx);
		visu.update(p);
	}

	return 0;
}
