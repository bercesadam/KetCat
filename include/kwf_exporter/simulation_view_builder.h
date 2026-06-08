#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atomic_units.h"
#include "laser/laser_pulse.h"
#include "simulation_data.h"

#include "local_space/local_space.h"
#include "local_space/density_matrix.h"
#include "global_space/output_probabilities.h"


namespace KetCat
{
    template <NeutralAtomTypeConfig Config, natural_t QubitCount>
    class SimulationViewBuilder
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;
        using OperationSpace = typename Manifold::SingleAtomOperationHilbertSpace;
        using FullHilbertSpace = typename Manifold::SingleAtomFullHilbertSpace;
        using DensityMatrixHelper = DensityMatrix<ConfigType::LevelCount, QubitCount>;

		/// @brief Keeping the name of the simulation steps for each qubit
		/// to be displayed for idle qubits in the visualization.
        /// (Simulation steps names are provided with the active gate name only.)
		std::array<std::string, QubitCount> m_lastGateName;

		/// @brief Reference to the NeutralAtomManifold, used to access the basis
        // states and projection methods for visualization.
        const Manifold& m_Manifold;

    public:
        SimulationViewBuilder(const Manifold& manifold)
        : m_Manifold(manifold)
        {
			std::fill(m_lastGateName.begin(), m_lastGateName.end(), "Idle");
        }

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
        
        SimulationView<FullHilbertSpace, QubitCount> build(
            const std::string simulationStepName,
            const real_t time,
            const StateVector<OperationSpace>& psiOp,
			const natural_t activeQubitIndex,
			const std::optional<natural_t> activeQubitIndex2,
            const LaserPulse& laser1,
            const LaserPulse& laser2) const
        {
            SimulationView<FullHilbertSpace> Data;

			// Simulation time in seconds (converted from atomic units)
            Data.m_time = Units::AtomicTimeToSeconds * time;

			// Get logical output probabilties from the full state vector
			Data.m_outputProbabilities = extractLogicalProbabilities<QubitCount, Config>(psiOp, true);

            for (natural_t i = 0; i < QubitCount; ++i)
            {
				Matrix<ConfigType::LevelCount> Rho = DensityMatrixHelper::reducedDensityMatrix(psiOp, i);

                // For each qubit, compute the purity of the reduced density matrix
                Data.m_qubitDatum[i].m_purity = DensityMatrixHelper::purity(Rho);

                // Project the local tile state back to the full Hilbert space for visualization
                Data.m_psi2D = m_Manifold.projectToFullHilbertSpace(Rho);

				// If this qubit is one of the active qubits, include the laser parameters in the data for visualization
				if (i == activeQubitIndex || (activeQubitIndex2 && i == *activeQubitIndex2))
				{
                    Data.m_laser1Wavelength = Units::wavelengthNmFromOmegaAu(laser1.m_omega);
                    Data.m_laser1Intensity = Units::intensityWcm2FromFieldAu(laser1.m_amplitude);
                    Data.m_laser2Wavelength = Units::wavelengthNmFromOmegaAu(laser2.m_omega);
                    Data.m_laser2Intensity = Units::intensityWcm2FromFieldAu(laser2.m_amplitude);

                    Data.m_title = compileCaption(simulationStepName, time, psiOp);
                    m_lastGateName[i] = simulationStepName;
				}
				// For non-active qubits, set laser parameters to zero
                else
                {
                    Data.m_laser1Wavelength = {};
                    Data.m_laser1Intensity  = {};
                    Data.m_laser2Wavelength = {};
                    Data.m_laser2Intensity  = {};
                    Data.m_title = compileCaption("Last gate: " + m_lastGateName[i], time, psiOp);
                }

            }

            return Data;
        }

    private:
        
    };

}
