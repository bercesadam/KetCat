#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief A classical quantum physics textbook example: a Gaussian wave packet in an infinite potential well
/// with a barrier in the middle (can be configured or removed to see free propagation).
/// This example demonstrates quantum tunneling through the potential barrier, which is generally the way
/// components of flash memories work in practice (electrons tunneling through oxide barriers).

int main()
{
	constexpr real_t BoxLength = 1.0;
	constexpr natural_t Steps = 96;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;

    constexpr KetCat::real_t x0 = 0.01; 
	constexpr KetCat::real_t sigma = 0.1;
	constexpr KetCat::real_t k0 = ConstexprMath::Pi * 10;
	constexpr auto GaussianPacket = FreeParticleGaussianWavePacket<HilbertSpace>()(x0, k0, sigma);

	constexpr PotentialBarrier potentialBarrier
	{
		0.45, 0.55, // potential wall in the middle
		3000
	};

	constexpr real_t TimeStep = 1E-4;
	constexpr KetCat::real_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<HilbertSpace::Dim>(mass, HilbertSpace::dx, potentialBarrier);

	OneDimensionalParticleBox<HilbertSpace> box(hamiltonian, GaussianPacket, TimeStep);
	
	
	Visu::VisuOscilloscope<HilbertSpace::Dim> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::NO,
		Visu::ShowPotential::YES
	);

	visu.setPotential(potentialBarrier, HilbertSpace::dx);

	while (true)
	{
		auto p = box.evolve();
		p.normalize();
		visu.update(p);
	}

	return 0;
}
