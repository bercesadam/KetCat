
#include <memory>
#include <functional>

#include "systems/time_master.h"
#include "systems/neutral_atom_manifold.h"
#include "laser/single_qbit_control.h"
#include "kwf_exporter/simulation_view_builder.h"
#include "kwf_exporter/kwf_exporter.h"

// Smoke test for the single-qubit control layer, applying a sequence of gates to a Cs atom and exporting the state evolution to a KWF file.

using namespace KetCat;

int main() 
{
    using namespace SpectroscopicLetters;

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

    TimeMaster::Clock().init(100);
    SingleQubitControl<Config> Controller;

    SimulationViewBuilder<Config> ViewBuilder(Manifold);
    StateVectorExporter<FullHilbertSpace> Exporter
    (
        "simulation.kwf",
        KetCat::ExportMode::RealImag
    );

	std::string SimuStep;
    natural_t FrameCounter = 0;
	natural_t SaveNthFrame = 2E6;

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
    /*
    SimuStep = "Gate: Pauli-X, Laser Protocol: STIRAP";
	std::cout << "Applying X gate..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::X, ConstexprMath::Pi }, Psi, ExporterCallback);

    SimuStep = "Gate: Pauli-X, Laser Protocol: Inverted STIRAP";
	std::cout << "Applying X gate (Inverted STIRAP)..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::X, ConstexprMath::Pi }, Psi,  ExporterCallback);
    */
    SimuStep = "Gate: Pauli-Y, Laser Protocol: STIRAP";
    std::cout << "Applying Y gate..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::Y, ConstexprMath::Pi }, Psi, ExporterCallback);
    // wasFractionalStirap false
    // !even

    // even
	// wasFractionalStirap false

    SimuStep = "Gate: Rx(Pi/2), Laser Protocol: Inverted STIRAP";
    std::cout << "Applying Rx(Pi/2)  gate..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::X, ConstexprMath::Pi / 2 }, Psi, ExporterCallback);
    // even
    // wasFractionalStirap false

    // even
    // wasFractionalStirap false


    std::cout << "Applying Z gate (virtual)..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::Z, ConstexprMath::Pi }, Psi, ExporterCallback);

    SimuStep = "Gate: Hadamard, Laser Protocol: STIRAP + Virtual Z";
    std::cout << "Applying Ry(Pi/2) gate..." << std::endl;
    Controller.applyPulseCommand({ RotationAxis::Y, ConstexprMath::Pi / 2 }, Psi, ExporterCallback);
}
