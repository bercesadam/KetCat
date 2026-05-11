#include "visu/matplot_based/visu_oscilloscope_3d.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief Sanity check: equal superposition of the first two eigenstates of a zero potential well

int main()
{
	constexpr real_t BoxLength = 1.0;
	constexpr natural_t Steps = 1024;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;

	constexpr auto Psi1 = EigenState<HilbertSpace>()(1).m_Psi;
	constexpr auto Psi2 = EigenState<HilbertSpace>()(2).m_Psi;
	complex_t alpha = complex_t::fromReal(1.0 / std::sqrt(2.0));
	auto Psi = Psi1.superpose(Psi2, alpha, alpha); // Superposition of the first two eigenstates with equal weights

	real_t Time = 0.0;
	constexpr real_t TimeStep = 3E-3;
	constexpr KetCat::real_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<HilbertSpace::Dim>(mass, HilbertSpace::dx, ZeroPotential);

	OneDimensionalParticleBox<HilbertSpace> box(hamiltonian, Psi, TimeStep);
	
	Visu::VisuOscilloscope3D visu;
	visu.setLimits(-3, 3, -3, 3);

	while (true)
	{
		Time += TimeStep;
		auto p = box.evolve();
		p.normalize();
		std::string Title = std::format("|ψ⟩ = 1/√2 (|1⟩ + |2⟩), t = {:.3f} a.u.", Time);
		visu.update(p, Title);
	}

	return 0;
}
