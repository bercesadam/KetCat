#include <iostream>
#include <iomanip>
#include <memory>

#include "systems/neutral_atom_manifold.h"
#include "laser/laser_drive_protocol.h"
#include "kwf_exporter/kwf_exporter.h"


using namespace KetCat;


int main() 
{
    using namespace SpectroscopicLetters;

    NeutralAtomTypeConfig
    <
        Element::Cs, 256, 75.0,
        0, 1, 4,
        QuantumNumber<6, s>,
        QuantumNumber<6, p>,
        QuantumNumber<7, s>,
        QuantumNumber<7, p>,
        QuantumNumber<40, s>,
        QuantumNumber<40, p>,
        QuantumNumber<39, s>,
        QuantumNumber<39, p>,
        QuantumNumber<5, d>,
        QuantumNumber<6, d>
    > Config;

	NeutralAtomManifold<Config> Manifold;
    using HilbertSpace = typename decltype(Manifold)::SingleAtomFullHilbertSpace;

    StateVectorExporter<HilbertSpace> Exporter
    (
        "simulation.kwf",
        KetCat::ExportMode::RealImag
    );

    auto ExporterCallback = [&](real_t time, const auto& currentPsi, const SiLaserPulse& laser1, const SiLaserPulse& laser2)
    {

        std::ostringstream Title;
        Title << std::fixed << std::setprecision(2);
        Title << "Cesium (Z=55)|";
        //Title << "Laser: " << wavelength_nm << " nm; " << requiredIntensity_Wcm2 << " W/cm²|";
        Title << "Time: " << time << " a.u.|";
        Title << "Populations: ";
        Title << "6s: " << currentPsi[0].normSquared() * 100.0 << "% ";
        Title << "6p: " << currentPsi[1].normSquared() * 100.0 << "% ";
        Title << "7s: " << currentPsi[2].normSquared() * 100.0 << "% ";
        Title << "7p: " << currentPsi[3].normSquared() * 100.0 << "% ";

        auto Psi2D = Manifold.projectToFullHilbertSpace(currentPsi);
        Exporter.writeTimestep(time, Psi2D, Title.str(), currentPsi[0], currentPsi[2], laser1, laser2);

        for (natural_t i = 0; i < Config.LevelCount; ++i)
        {
            std::cout << "Probability of basis state" << i << ": " << currentPsi[i].normSquared() * 100.0 << "%" << std::endl;

        }
        std::cout << "------------------------" << std::endl;
    };

    auto Psi = Manifold.getOperationSeed();
    LaserDrivePulse<Config> STIRAP;
    STIRAP(Psi, Config.Logical0Level, { RotationAxis::X, ConstexprMath::Pi / 2 }, ExporterCallback);
    STIRAP(Psi, Config.Logical0Level, { RotationAxis::Z, ConstexprMath::Pi / 2 }, ExporterCallback);
    STIRAP(Psi, Config.Logical0Level, { RotationAxis::Y, ConstexprMath::Pi / 2 }, ExporterCallback);
}
