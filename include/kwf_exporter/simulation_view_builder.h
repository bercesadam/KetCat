#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atomic_units.h"
#include "laser/laser.h"
#include "simulation_data.h"


namespace KetCat
{
    template <NeutralAtomTypeConfig Config>
    class SimulationViewBuilder
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;

    public:
        using Manifold = NeutralAtomManifold<Config>;
        using OperationSpace = typename Manifold::SingleAtomOperationHilbertSpace;
        using FullHilbertSpace = typename Manifold::SingleAtomFullHilbertSpace;

        SimulationViewBuilder(const Manifold& manifold)
        : m_Manifold(manifold)
        {}

        std::string compileCaption(const std::string simulationStepName,
            real_t time, const StateVector<OperationSpace>& psiOp) const
        {
            std::ostringstream Title;
            Title << std::fixed << std::setprecision(2);
            Title << "Atom: " << getElementName(ConfigType::ChemicalElement);
            Title << " (Z=" << AtomicNumber<ConfigType::ChemicalElement>::value << ")|";
            Title << simulationStepName << "|"; 
            Title << "Time: " << Units::AtomicTimeToSeconds * time * 1E9 << " ns|";

            Title << "Populations: {";

            [&] <size_t... IndexSequence>(std::index_sequence<IndexSequence...>) {
                ((Title << (IndexSequence > 0 ? ", " : "")
                    << quantumNumberToString<std::tuple_element_t<IndexSequence, typename ConfigType::QuantumNumbers>>()
                    << ": " << psiOp[IndexSequence].normSquared() * 100.0 << "%"), ...);
            }(std::make_index_sequence<ConfigType::LevelCount>{});

            Title << "}";

			return Title.str();
        }
        
        SimulationView<FullHilbertSpace> build(
            const std::string simulationStepName,
            const real_t time,
            const StateVector<OperationSpace>& psiOp,
            const LaserPulse& laser1,
            const LaserPulse& laser2) const
        {
            SimulationView<FullHilbertSpace> Data;
            Data.m_time = Units::AtomicTimeToSeconds * time;

            Data.m_psi2D = m_Manifold.projectToFullHilbertSpace(psiOp);

            Data.m_alpha = psiOp[ConfigType::Logical0Level];
            Data.m_beta  = psiOp[ConfigType::Logical1Level];

            Data.m_laser1Wavelength = Units::wavelengthNmFromOmegaAu(laser1.m_omega);
            Data.m_laser1Intensity = Units::intensityWcm2FromFieldAu(laser1.m_amplitude);
            Data.m_laser2Wavelength = Units::wavelengthNmFromOmegaAu(laser2.m_omega);
            Data.m_laser2Intensity = Units::intensityWcm2FromFieldAu(laser2.m_amplitude);

            Data.m_title = compileCaption(simulationStepName, time, psiOp);
            return Data;
        }

    private:
        const Manifold& m_Manifold;
    };

}
