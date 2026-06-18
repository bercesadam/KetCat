#pragma once
#include "core_types.h"

#include "quantum_processor/time_master.h"

#include "local_space/neutral_atom_manifold.h"
#include "global_space/subspace_operations.h"

#include "kwf_exporter/simulation_view_builder.h"
#include "kwf_exporter/kwf_exporter.h"


namespace KetCat
{
    /// @brief Constant used to force an export step regardless of the frame decimation counter.
    /// Makes exportStep() calls look better, but I found defining an enum class a bit of an overkill.
    static constexpr bool KEYFRAME = true;

    static constexpr bool VERBOSE = true;

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
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;
        using VisuHilbertSpace = typename NeutralAtomManifold<Config>::SingleAtomFullHilbertSpace;
        using FullHilbertSpace = typename GlobalStateManager::FullHilbertSpace;

        /// @brief Transformer to map numerical simulation data to visualizable structures.
        SimulationViewBuilder<Config, QubitCount> m_ViewBuilder;

        /// @brief Binary stream handler for simulation data export
        StateVectorExporter<VisuHilbertSpace, QubitCount> m_Exporter;

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
              m_Exporter(fileName, QubitCount), 
              m_SaveNthFrame(saveNthFrame)
        {
        }

        /// @brief Update the metadata label for the current sequence of frames.
        void appendSimulationStepName(const std::string& stepName)
        {
            m_SimulationStepName += stepName;
        }

        //  @brief Empty the metadata label 
        void resetSimulationStepName()
        {
            m_SimulationStepName.clear();
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
        void exportStep(const StateVector<FullHilbertSpace, QuantumPicture::Schrodinger>& psi,
            const qbit_list_t<2> targets, const LaserPulse& laser1, const LaserPulse& laser2,
            const bool isKeyFrame = false)
        {
            if (m_FrameCounter % m_SaveNthFrame == 0 || isKeyFrame)
            {
                auto SimulationView =
                    m_ViewBuilder.build(
                        m_SimulationStepName,
                        TimeMaster::Clock().getGlobalTime(),
                        psi, targets, laser1, laser2);

                m_Exporter.writeTimestep(SimulationView);

                if (isKeyFrame || VERBOSE)
                {
                    std::cout << "State vector: " << psi << std::endl;
                    for (natural_t q = 0; q < QubitCount; ++q)
                    {
						std::cout << SimulationView.m_qubitDatum[q].m_title << std::endl;
						std::cout << "Purity of qubit " << q << ": " << SimulationView.m_qubitDatum[q].m_purity << std::endl;
                    }

                    std::cout << "------------------------" << std::endl << std::endl;
                }
            }

            m_FrameCounter++;
        }

    };
}