
#include <memory>
#include <functional>

#include "quantum_computer/time_master.h"
#include "quantum_bit/neutral_atom_manifold.h"
#include "laser/single_qbit_control.h"
#include "kwf_exporter/simulation_view_builder.h"
#include "kwf_exporter/kwf_exporter.h"

// Smoke test for the single-qubit control layer, applying a sequence of gates to a Cs atom and exporting the state evolution to a KWF file.

using namespace KetCat;

int main() 
{
    using namespace SpectroscopicLetters;

	real_t TimeStep = 50;

    NeutralAtomTypeConfig
    <
        Element::Cs,

        256, /* Spatial discretization steps count */
        100.0, /* Spatial extent in a.u. */

        0, /* Index of the logical level 0 */
        2, /* Index of the logical level 1*/
		4, /* Index of the Rydberg level */

        QuantumNumber<6, s>,
        QuantumNumber<6, p>,
        QuantumNumber<7, s>,
        QuantumNumber<7, p>,
        QuantumNumber<40, s>,
        QuantumNumber<40, p>
    > Config;

	NeutralAtomManifold<Config> Manifold;
    using FullHilbertSpace = typename decltype(Manifold)::SingleAtomFullHilbertSpace;
	using OperationHilbertSpace = typename decltype(Manifold)::SingleAtomOperationHilbertSpace;

    TimeMaster::Clock().init(TimeStep);
    SingleQubitControl<Config> Controller;

    SimulationViewBuilder<Config> ViewBuilder(Manifold);
    StateVectorExporter<FullHilbertSpace> Exporter
    (
        "simulation.kwf",
        KetCat::ExportMode::RealImag
    );

	std::string SimuStep;
    natural_t FrameCounter = 0;
	natural_t SaveNthFrame = 1E6;

	decltype(Controller)::CallbackType ExporterCallback =
        [&](const StateVector<OperationHilbertSpace>& currentPsi, const LaserPulse& laser1, const LaserPulse& laser2, const bool isKey)
        {
            if (FrameCounter % SaveNthFrame == 0 || isKey)
            {

                auto SimulationView = ViewBuilder.build(SimuStep, TimeMaster::Clock().getGlobalTime(), currentPsi, laser1, laser2);
                Exporter.writeTimestep(SimulationView);

                for (natural_t i = 0; i < Config.LevelCount; ++i)
                {
                    std::cout << "Probability of basis state " << i << ": " << currentPsi[i].normSquared() * 100.0 << "%" << std::endl;
                }
                std::cout << "------------------------" << std::endl;
            }
            FrameCounter++;
        };

    StateVector<OperationHilbertSpace> Psi = Manifold.getOperationSeed();

    std::cout << "Applying X gate..." << std::endl;
    SimuStep = "Gate: Pauli-X, Laser Protocol: STIRAP";
    Controller.applyPulseCommand({ RotationAxis::X, ConstexprMath::Pi }, Psi, ExporterCallback);

    std::cout << "Applying Y gate..." << std::endl;
    SimuStep = "Gate: Pauli-Y, Laser Protocol: Inverted STIRAP";
    Controller.applyPulseCommand({ RotationAxis::Y, ConstexprMath::Pi }, Psi, ExporterCallback);

    std::cout << "Applying H gate..." << std::endl;
    SimuStep = "Gate: Hadamard [Virtual Z(Pi) + Ry(Pi/2)], Laser Protocol: Fractional STIRAP";
    Controller.applyPulseCommand({ RotationAxis::Z, ConstexprMath::Pi }, Psi, ExporterCallback);
    Controller.applyPulseCommand({ RotationAxis::Y, ConstexprMath::Pi / 2 }, Psi, ExporterCallback);
}
