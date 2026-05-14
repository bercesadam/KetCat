#pragma once
#include "core_types.h"

#include "quantum_processor/time_master.h"

#include "operation_space/neutral_atom_manifold.h"
#include "operation_space/utils/subspace_operations.h"

#include "kwf_exporter/simulation_view_builder.h"
#include "kwf_exporter/kwf_exporter.h"


namespace KetCat
{
    /// @brief Constant used to force an export step regardless of the frame decimation counter.
    /// Makes exportStep() calls look better, but I found defining an enum class a bit of an overkill.
    constexpr bool KEYFRAME = true;

    /// @brief Diagnostic and telemetry layer for recording simulation state.
    /// 
    /// @details
    ///    The SimulationObserver captures the high-dimensional quantum state and 
    ///    associated laser control parameters, serializing them into a binary format (.kwf).
    ///
    ///    Key Responsibilities:
    ///      • Decimation: Filters high-frequency TDSE steps to a manageable framerate.
    ///      • State Mapping: Converts operational subspace vectors to full-manifold (2D) wavefunctions.
    ///      • Persistence: Interfaces with the KWF exporter for binary storage.
    /// 
    /// @tparam QubitCount Number of qubits in the register.
    /// @tparam Config Atomic configuration defining the mapping and energy levels.
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class SimulationObserver
    {
        using FullHilbertSpace = typename NeutralAtomManifold<Config>::SingleAtomFullHilbertSpace;
        using OperationHilbertSpace = typename NeutralAtomManifold<Config>::SingleAtomOperationHilbertSpace;

        /// @brief Transformer to map numerical simulation data to visualizable structures.
        SimulationViewBuilder<Config> m_ViewBuilder;

        /// @brief Binary stream handler for simulation data persistence.
        StateVectorExporter<FullHilbertSpace> m_Exporter;

        /// @brief Label for the current logical operation being recorded.
        std::string m_SimulationStepName;

        /// @brief Frequency of data capture (records every N-th frame).
        natural_t m_SaveNthFrame;

        /// @brief Internal counter for frame decimation logic.
        natural_t m_FrameCounter = 0;

    public:
        /// @brief Initialize the observer with a target file and capture frequency.
        ///
        /// @param manifold Reference to the physical atomic manifold.
        /// @param fileName Path to the output .kwf file.
        /// @param saveNthFrame Frame decimation interval.
        SimulationObserver(const NeutralAtomManifold<Config>& manifold,
            const std::string fileName, const natural_t saveNthFrame)
            : m_ViewBuilder(manifold),
            m_Exporter(fileName, ExportMode::RealImag),
            m_SaveNthFrame(saveNthFrame)
        {
            std::cout << "SimulationObserver initialized with file: " << fileName
				<< " and save frequency: " << saveNthFrame << std::endl;
        }

        /// @brief Update the metadata label for the current sequence of frames.
        void setSimulationStepName(const std::string& stepName)
        {
            m_SimulationStepName = stepName;
        }

        /// @brief Sample the current quantum and laser state for export.
        ///
        /// @param psi The current single-qubit operational wavefunction.
        /// @param laser1 Current state of the pump laser.
        /// @param laser2 Current state of the Stokes laser.
        /// @param isKeyFrame If true, ignores m_SaveNthFrame and forces a write.
        ///
        /// @details
        ///    Calculates the instantaneous basis state probabilities and passes 
        ///    the formatted view to the binary exporter if the capture criteria are met.
        void exportStep(const StateVector<OperationHilbertSpace>& psi,
            const LaserPulse& laser1, const LaserPulse& laser2,
            const bool isKeyFrame = false)
        {
            if (m_FrameCounter % m_SaveNthFrame == 0 || isKeyFrame)
            {
                auto SimulationView =
                    m_ViewBuilder.build(
                        m_SimulationStepName, TimeMaster::Clock().getGlobalTime(),
                        psi, laser1, laser2);

                m_Exporter.writeTimestep(SimulationView);

                std::cout << "Exported frame " << m_FrameCounter << ": " << m_SimulationStepName
					<< " at time " << SimulationView.m_time * 1E9 << " ns" << std::endl;

                // Diagnostic terminal output, only for debugging, to be prettified or removed
                for (natural_t i = 0; i < decltype(Config)::LevelCount; ++i)
                {
                    std::cout << "Probability of basis state " << i << ": " << psi[i].normSquared() * 100.0 << "%" << std::endl;
                }
                std::cout << "------------------------" << std::endl;
            }

            m_FrameCounter++;
        }

    };
}