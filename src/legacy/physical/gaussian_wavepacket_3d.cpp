#include <format>
#include "visu/matplot_based/visu_oscilloscope_3d.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;
using namespace KetCat::Visu;

///@brief A classical quantum physics textbook example: a Gaussian wave packet in an infinite potential well
/// with a barrier in the middle (can be configured or removed to see free propagation).
/// This example demonstrates quantum tunneling through the potential barrier, which is generally the way
/// components of flash memories work in practice (electrons tunneling through oxide barriers).

int main()
{
	constexpr real_t BoxLength = 3.0;
	constexpr natural_t Steps = 2048;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;

    constexpr KetCat::real_t x0 = 1.5; 
	constexpr KetCat::real_t sigma = 0.3;
	constexpr KetCat::real_t k0 = ConstexprMath::Pi * 5.0;
	constexpr auto GaussianPacket = FreeParticleGaussianWavePacket<HilbertSpace>()(x0, k0, sigma);

	real_t Time = 0.0;
	constexpr real_t TimeStep = 1E-3;
	constexpr KetCat::real_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<HilbertSpace::Dim>(mass, HilbertSpace::dx, ZeroPotential);

	OneDimensionalParticleBox<HilbertSpace> box(hamiltonian, GaussianPacket, TimeStep);
	
	
	Visu::VisuOscilloscope3D visu;
	visu.setLimits(-5, 5, -5, 5);
	unsigned Frame = 0;

	while (true)
	{
		Time += TimeStep;
		auto p = box.evolve();
		p.normalize();

		std::string Title = std::format("Gaussian wave packet, t = {:.3f} a.u.", Time);
		visu.update(p, Title);
	}

	return 0;
}
