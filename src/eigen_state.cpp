#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

///@brief he selected eigenstate forms a standing wave in space, meaning the probability density remains stationary.
/// However, the real and imaginary components of the wavefunction continuously rotate in the complex plane due to
/// the time-dependent phase factor. As a result, while the spatial profile appears fixed (standing wave),
/// the Re and Im parts oscillate in time, demonstrating the internal complex phase rotation of the quantum state.

int main()
{
	constexpr real_t BoxLength = 1.0;
	constexpr natural_t Steps = 136;
	using HilbertSpace = InfiniteHilbertSpace<1_D, Steps, BoxLength>;

	constexpr auto Psi = EigenState<HilbertSpace>()(5).m_Psi;

	constexpr real_t TimeStep = 1E-3;
	constexpr KetCat::real_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<HilbertSpace::Dim>(mass, HilbertSpace::dx, ZeroPotential);

	OneDimensionalParticleBox<HilbertSpace> box(hamiltonian, Psi, TimeStep);
	
	Visu::VisuOscilloscope<HilbertSpace::Dim> visu(
		Visu::UsePhaseEncoding::YES,
		Visu::ClearScreen::YES,
		Visu::ShowComplexParts::YES,
		Visu::ShowPotential::NO
	);

	while (true)
	{
		auto p = box.evolve();
		p.normalize();
		visu.update(p);
	}

	return 0;
}
