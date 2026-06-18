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
		using GlobalSpace = typename SubspaceHelper<ConfigType::LevelCount, QubitCount>::FullHilbertSpace;
        using DensityMatrixHelper = DensityMatrix<ConfigType::LevelCount, QubitCount>;

		/// @brief Keeping the name of the simulation steps for each qubit
		/// to be displayed for idle qubits in the visualization.
        /// (Simulation steps names are provided with the active gate name only.)
		std::array<std::string, QubitCount> m_lastGateName;

		/// @brief Reference to the NeutralAtomManifold, used to access the basis
        // states and projection methods for visualization.
        const Manifold& m_Manifold;

		/// @brief Helper class to extract logical probabilities from the full global state vector.
		LogicalProbabilityExtractor<QubitCount, Config> m_ProbabilityExtractor;

    public:
        SimulationViewBuilder(const Manifold& manifold)
        : m_Manifold(manifold)
        {
			std::fill(m_lastGateName.begin(), m_lastGateName.end(), "Idle");
        }

        std::string compileCaption(const std::string simulationStepName,
            real_t time, const StateVector<OperationSpace, QuantumPicture::Schrodinger>& psiLocal) const
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
                    << ": " << psiLocal[IndexSequence].normSquared() * 100.0 << "%"), ...);
            }(std::make_index_sequence<ConfigType::LevelCount>{});

            Title << "}";

			return Title.str();
        }
        
        SimulationView<FullHilbertSpace, QubitCount> build(
            const std::string simulationStepName,
            const real_t time,
            const StateVector<GlobalSpace, QuantumPicture::Schrodinger>& psiOp,
            const qbit_list_t<2> targets,
            const LaserPulse& laser1,
            const LaserPulse& laser2)
        {
            SimulationView<FullHilbertSpace, QubitCount> Data;

			// Simulation time in seconds (converted from atomic units)
            Data.m_time = Units::AtomicTimeToSeconds * time;

			// Get logical output probabilties from the full state vector
			Data.m_outputProbabilities =
                m_ProbabilityExtractor.extractLogicalProbabilities(psiOp, true);

			std::cout << "Probabilities: ";
			for (size_t i = 0; i < Data.m_outputProbabilities.size(); ++i)
			{
				std::cout << "P(" << std::bitset<QubitCount>(i) << ")=" << Data.m_outputProbabilities[i] * 100.0 << "% ";
			}
            std::cout << std::endl; 

            for (natural_t i = 0; i < QubitCount; ++i)
            {
                QubitData<FullHilbertSpace>& Qubit = Data.m_qubitDatum[i];

				Matrix<ConfigType::LevelCount> Rho = DensityMatrixHelper::reducedDensityMatrix(psiOp, i);

                // For each qubit, compute the purity of the reduced density matrix
                Qubit.m_purity = DensityMatrixHelper::purity(Rho);

				// Convert the reduced density matrix to a state vector in the local operation space for visualization
				StateVector<OperationSpace, QuantumPicture::Schrodinger> psiLocal = m_Manifold.getReducedStateFromDensityMatrix(Rho);

				// Extract the logical amplitudes (|0⟩ and |1⟩) for the qubit from the local state vector
                Qubit.m_alpha = psiLocal[ConfigType::Logical0Level];
				Qubit.m_beta  = psiLocal[ConfigType::Logical1Level];

                // Project the local tile state back to the full Hilbert space for visualization
                Qubit.m_psi2D = m_Manifold.projectToFullHilbertSpace(psiLocal);

				// If this qubit is one of the active qubits, include the laser parameters and update title
				if (i == targets[0] || i == targets[1])
				{
                    Qubit.m_laser1Wavelength = Units::wavelengthNmFromOmegaAu(laser1.m_omega);
                    Qubit.m_laser1Intensity = Units::intensityWcm2FromFieldAu(laser1.m_amplitude);
                    Qubit.m_laser2Wavelength = Units::wavelengthNmFromOmegaAu(laser2.m_omega);
                    Qubit.m_laser2Intensity = Units::intensityWcm2FromFieldAu(laser2.m_amplitude);

                    Qubit.m_title =
                        compileCaption("Qubit " + std::to_string(i) + " " + simulationStepName, time, psiLocal);
                    m_lastGateName[i] = simulationStepName;
				}
				// For non-active qubits, set laser parameters to zero and print the last active gate name
                else
                {
                    Qubit.m_laser1Wavelength = {};
                    Qubit.m_laser1Intensity  = {};
                    Qubit.m_laser2Wavelength = {};
                    Qubit.m_laser2Intensity  = {};
                    Qubit.m_title =
                        compileCaption("Qubit " + std::to_string(i) + " is idle, last gate: " + m_lastGateName[i], time, psiLocal);
                }

            }

            return Data;
        }

    private:
        
    };

}
