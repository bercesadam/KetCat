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
              m_Exporter(fileName), 
              m_SaveNthFrame(saveNthFrame)
        {
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
        void exportStep(const StateVector<FullHilbertSpace>& psi,
            const LaserPulse& laser1, const LaserPulse& laser2,
            const bool isKeyFrame = false)
        {
            if (m_FrameCounter % m_SaveNthFrame == 0 || isKeyFrame)
            {
                // Print full state vector
                if (isKeyFrame)
                {
                    for (natural_t i = 0; i < GlobalStateManager::FullDim; ++i)
                    {
                        std::cout << "Global basis state " << i << ": Re: " << psi[i].re << "\tIm: " << psi[i].im << std::endl;
                    }
                    std::cout << "------------------------" << std::endl << std::endl;
                }

                for (natural_t q = 0; q < QubitCount; ++q)
                {
                    auto qubitLocalState = GlobalStateManager::extractLocalState(psi, q);

                    auto SimulationView =
                        m_ViewBuilder.build(
                            "Qubit" + std::to_string(q) + ": " + m_SimulationStepName,
                            TimeMaster::Clock().getGlobalTime(),
                            psi, laser1, laser2);

                    m_Exporter.writeTimestep(SimulationView);

                    if (isKeyFrame)
                    {
                        std::cout << "Purity of qubit " << q << ": " << qubitLocalState.purityValue << std::endl;
                        //std::cout << "Laser intensities in au: Pump = " << laser1.m_amplitude << ", Stokes = " << laser2.m_amplitude << std::endl;  

                        for (natural_t i = 0; i < decltype(Config)::LevelCount; ++i)
                        {
                            std::cout << "Probability of basis state " << i << ": " << qubitLocalState.pureStateVector[i].normSquared() * 100.0 << "%\t";
                            std::cout << "Re: " << qubitLocalState.pureStateVector[i].re << "\tIm: " << qubitLocalState.pureStateVector[i].im << std::endl;
                        }


                        std::cout << "Exported timeframe " << m_FrameCounter << ": " << m_SimulationStepName
                            << " at time " << TimeMaster::Clock().getGlobalTime() * Units::AtomicTimeToSeconds * 1E9 << " ns" << std::endl;
                        std::cout << "------------------------" << std::endl << std::endl;
                    }
                }

                m_FrameCounter++;
            }
        }

    };
}