
#include <memory>
#include <functional>

#include "systems/neutral_atom_manifold.h"
#include "laser/single_qbit_protocol.h"
#include "kwf_exporter/simulation_view_builder.h"
#include "kwf_exporter/kwf_exporter.h"


using namespace KetCat;


int main() 
{
    using namespace SpectroscopicLetters;

    NeutralAtomTypeConfig
    <
        Element::Cs, 256, 75.0,
        0, 2, 4,
        QuantumNumber<6, s>,
        QuantumNumber<6, p>,
        QuantumNumber<7, s>,
        QuantumNumber<7, p>,
        QuantumNumber<40, s>/*,
        QuantumNumber<40, p>,
        QuantumNumber<39, s>,
        QuantumNumber<39, p>,
        QuantumNumber<5, d>,
        QuantumNumber<100, f>*/
    > Config;

	NeutralAtomManifold<Config> Manifold;
    using HilbertSpace = typename decltype(Manifold)::SingleAtomFullHilbertSpace;
	using OperationHilbertSpace = typename decltype(Manifold)::SingleAtomOperationHilbertSpace;

    SimulationViewBuilder<Config> ViewBuilder(Manifold);
    StateVectorExporter<HilbertSpace> Exporter
    (
        "simulation.kwf",
        KetCat::ExportMode::RealImag
    );

    natural_t FrameCounter = 0;
    std::function<void(real_t, const StateVector<OperationHilbertSpace>&, const LaserPulse&, const LaserPulse&)> ExporterCallback =
        [&](real_t time, const StateVector<OperationHilbertSpace>& currentPsi, const LaserPulse& laser1, const LaserPulse& laser2)
        {
            if (FrameCounter % 500000 == 0)
            {
                auto SimulationView = ViewBuilder.build("Gate: Pauli-X", time, currentPsi, laser1, laser2);
                Exporter.writeTimestep(SimulationView);

                for (natural_t i = 0; i < Config.LevelCount; ++i)
                {
                    std::cout << "Probability of basis state" << i << ": " << currentPsi[i].normSquared() * 100.0 << "%" << std::endl;
                }
                std::cout << "------------------------" << std::endl;
            }

            FrameCounter++;
        };

    StateVector<OperationHilbertSpace> Psi = Manifold.getOperationSeed();
    SingleQubitProtocol<Config> Protocol;
	Protocol.applyPulseCommand({ RotationAxis::X, ConstexprMath::Pi }, Psi, ExporterCallback);
}
