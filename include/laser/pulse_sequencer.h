#pragma once
#include <optional>
#include "two_photon_laser.h"
#include "quantum_processor/time_master.h"
#include "compiler/physical_instruction.h"
#include "operation_space/neutral_atom_manifold.h"

namespace KetCat
{
    namespace detail
    {
        /// @brief Internal state tracking for individual atom control.
        /// 
        /// @details 
        ///    Maintains the phase of the rotating frame and tracking for 
        ///    fractional STIRAP cycles to ensure adiabatic continuity.
        struct AtomControlData
        {
            /// @brief Current phase of the rotating frame φ_frame.
            real_t m_framePhase = 0.0;

            /// @brief Current position in the continuous STIRAP cycle θ_accum.
            /// @details Tracks the total mixing angle in the range [0, 2π).
            real_t m_stirapCycleTheta = 0.0;
              
			/// @brief Flag indicating if the atom is currently in the Rydberg state.
			/// @details As Rydberg excitations are always complete STIRAP cycles
			/// and also it's independent from the logical state, we  track this separately
            /// to don't mix with the concept behind m_stirapCycleTheta, but still
			/// allowing the control logic to know if STIRAP or Inverted STIRAP should be used
			bool m_isRydbergExcited = false;  

            bool m_isInitialized = false;

            /// @brief Time required for a complete adiabatic π-pulse.
            real_t m_stirapFullTransferTime = 0.0;

            /// @brief Current temporal offset within the pulse sequence.
            real_t m_stirapCycleTime = 0.0;

            /// @brief Select STIRAP protocol based on accumulated rotation angle.
            ///
            /// @param targetTheta
            ///    Requested logical rotation angle in radians (Δθ).
            ///
            /// @return
            ///    Selected two-photon transfer protocol.
            ///
            /// @details
            ///    Protocol selection alternates to maintain adiabaticity:
            ///
            ///      • STIRAP:          |0⟩ → |2⟩ when θ_accum ∈ [0, π)
            ///      • Inverted STIRAP: |2⟩ → |0⟩ when θ_accum ∈ [π, 2π)
            ///
            ///    The selection is governed by the state of the mixing angle:
            ///      tan(θ/2) = Ωp(t) / Ωs(t)
            TwoPhotonProtocol handleStirapTheta(real_t targetTheta)
            {
                TwoPhotonProtocol Result;

                if (m_stirapCycleTheta >= 0.0 &&
                    m_stirapCycleTheta < ConstexprMath::Pi)
                {
                    Result = TwoPhotonProtocol::STIRAP;
                }
                else if (m_stirapCycleTheta >= ConstexprMath::Pi - 1E-7)
                {
                    Result = TwoPhotonProtocol::InvertedSTIRAP;
                }

                m_stirapCycleTheta += targetTheta;

                /// Wrap accumulated angle: θ_accum ← θ_accum mod 2π
                if (ConstexprMath::floatNear(
                    m_stirapCycleTheta,
                    2 * ConstexprMath::Pi))
                {
                    m_stirapCycleTheta = 0.0;
                }

                return Result;
            }

            /// @brief Advance the internal pulse clock.
            void updateCycleTime(const real_t sequenceTime)
            {
                m_stirapCycleTime += sequenceTime;

                if (m_stirapCycleTime >= m_stirapFullTransferTime)
                {
                    m_stirapCycleTime = 0.0;
                }
            }

            void setFullTransferTime(const real_t fullTransferTime)
            {
                if (!m_isInitialized)
                {
                    m_stirapFullTransferTime = fullTransferTime;
                    m_isInitialized = true;
                }
            }
        };
    }

    /// @brief Single-qubit quantum control layer for Raman/STIRAP-driven operations.
    ///
    /// @details
    ///    Translates abstract PulseCommands into laser-driven dynamics.
    ///
    ///    Effective Raman Hamiltonian in the RWA:
    ///      H_eff = (ℏ/2) [ Ω_eff exp(iφ) |0⟩⟨2| + h.c. ]
    ///
    ///    Where:
    ///      Ω_eff ≈ (Ωp Ωs) / (2Δ)
    ///      φ_total = φ_axis - φ_frame
    ///
    /// @tparam Config Neutral atom level/dipole configuration.
    /// @tparam AtomCount Number of atoms in the register.
    template <NeutralAtomTypeConfig Config, natural_t AtomCount>
    class LaserPulseSequencer
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;

        const std::array<real_t, ConfigType::LevelCount> m_energies =
            Manifold::getHartreeEnergies();

        const square_matrix_t<ConfigType::LevelCount> m_dipoleMatrix =
            Manifold::getDipoleMatrix();

        real_t m_peakRabiHz = 50e6;
        real_t m_commonDetuningHz = 0;

        std::array<detail::AtomControlData, AtomCount> m_AtomControlState;


        /// @brief Construct physical laser configuration for a pulse command.
        ///
        /// @details
        ///    Calculates the required laser frequencies and phases.
        ///    ωp = (E₁ - E₀) + Δ
        ///    ωs = (E₂ - E₁) - Δ
        ///
        ///    The rotating frame update ensures phase coherence:
        ///      φ_pump = φ_instr - φ_frame
        TwoPhotonConfig prepareLaserConfig(detail::AtomControlData& atomControl, const PhysicalInstruction& instruction)
        {
            TwoPhotonConfig LaserConfig;

            LaserConfig.m_GroundLevelIndex = (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
				? ConfigType::Logical1Level
				: ConfigType::Logical0Level;

            LaserConfig.m_Level1Energy = m_energies[LaserConfig.m_GroundLevelIndex];
            LaserConfig.m_Level2Energy = m_energies[LaserConfig.m_GroundLevelIndex + 1];
            LaserConfig.m_Level3Energy = m_energies[LaserConfig.m_GroundLevelIndex + 2];

            LaserConfig.m_Mu12 = m_dipoleMatrix[LaserConfig.m_GroundLevelIndex][LaserConfig.m_GroundLevelIndex + 1].re;
            LaserConfig.m_Mu23 = m_dipoleMatrix[LaserConfig.m_GroundLevelIndex + 1][LaserConfig.m_GroundLevelIndex + 2].re;

			std::cout << "Hartree energy differences: " << std::endl;
			std::cout << "E1 - E0: " << LaserConfig.m_Level2Energy - LaserConfig.m_Level1Energy << " a.u." << std::endl;
			std::cout << "E2 - E1: " << LaserConfig.m_Level3Energy - LaserConfig.m_Level2Energy << " a.u." << std::endl;
			std::cout << "Dipole μ12: " << LaserConfig.m_Mu12 << " a.u." << std::endl;
			std::cout << "Dipole μ23: " << LaserConfig.m_Mu23 << " a.u." << std::endl;

            LaserConfig.m_peakRabiFrequency = Units::omegaAuFromHz(m_peakRabiHz);
            LaserConfig.m_commonDetuning = Units::omegaAuFromHz(m_commonDetuningHz);

            // Apply rotating frame correction: φ_eff = φ_laser - φ_frame
            LaserConfig.m_pumpPhase = instruction.m_phase - atomControl.m_framePhase;
            
            // For performing Rydberg blockades, which is independent from the control state
            // betweem the energy levels corresponding to the logical states we  toggle between
			// STIRAP and Inverted STIRAP protocols based on the excitation flag stored in atom control data.
            if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
            {
                if (!atomControl.m_isRydbergExcited)
                {
					LaserConfig.m_protocol = TwoPhotonProtocol::STIRAP;
                    atomControl.m_isRydbergExcited = true;
                }
                else
                {
					LaserConfig.m_protocol = TwoPhotonProtocol::InvertedSTIRAP;
                    atomControl.m_isRydbergExcited = false;
                }
            }
			// For performing one qubit logical rotations (X/Y), we toggle between STIRAP and Inverted STIRAP protocols
            // based on the accumulated rotation angle stored in atom control data.
            else
            {
                LaserConfig.m_protocol = atomControl.handleStirapTheta(instruction.m_theta);
            }

			// In case in inverted STIRAP sequences we need to adjust the target theta to be the
            // complementary angle to ensure the correct rotation direction
            LaserConfig.m_targetTheta = instruction.m_theta;
            if (LaserConfig.m_protocol == TwoPhotonProtocol::InvertedSTIRAP &&
                !(ConstexprMath::floatNear(ConstexprMath::Pi, LaserConfig.m_targetTheta)))
            {
                LaserConfig.m_targetTheta = ConstexprMath::Pi - instruction.m_theta;
            }

            return LaserConfig;
        }

    public:
        constexpr void initializeAtomAsLogical1(natural_t index)
        {
            m_AtomControlState[index].m_stirapCycleTheta = ConstexprMath::Pi;
        }

        /// @brief Generate the time-dependent envelope for a physical instruction.
        ///
        /// @param instruction The logical gate/pulse to be implemented.
        ///
        /// @return Optional envelope evaluator (empty for virtual Z).
        ///
        /// @details
        ///    Virtual Z rotations update the frame without physical pulses:
        ///      φ_frame ← φ_frame + θ
        ///
        ///    X/Y rotations use the TwoPhotonPulseBuilder to solve:
        ///      θ = ∫ Ω_eff(t) dt
        std::optional<TwoPhotonLaserEnvelope> calculateLaserEnvelope(const PhysicalInstruction& instruction)
        {
            const natural_t targetIndex = instruction.m_targets[0];
            if (targetIndex >= AtomCount)
            {
                return std::nullopt;
            }

            detail::AtomControlData& AtomControl = m_AtomControlState[targetIndex];

            // Handle Bloch-Z rotations (Virtual Z)
            if (instruction.m_type == PhysicalInstructionType::VirtualZ)
            {
                AtomControl.m_framePhase += instruction.m_theta;
                return std::nullopt;
            }

            // Generate physical parameters via Builder
            TwoPhotonPulseBuilder LaserBuilder(prepareLaserConfig(AtomControl, instruction));
            TwoPhotonLaserEnvelope LaserEnvelope = LaserBuilder.build();


            // For X/Y rotations, the instruction phase defines the axis
            // Update frame phase if relative phase tracking is required
            AtomControl.m_framePhase += instruction.m_phase;

            // Synchronize timing for adiabatic continuity
            // t_start allows resuming pulses from previous fractional steps
            LaserEnvelope.setStartTime(AtomControl.m_stirapCycleTime);

            AtomControl.setFullTransferTime(LaserEnvelope.getFullTransferTime());
            AtomControl.updateCycleTime(LaserEnvelope.getTransitionTimeLimit());

            return LaserEnvelope;
        }

        /// @brief Set peak Rabi frequency Ω_max / (2π).
        void setPeakRabiHz(real_t hz)
        {
            m_peakRabiHz = hz;
        }

        /// @brief Set common intermediate detuning Δ / (2π).
        void setDetuningHz(real_t hz)
        {
            m_commonDetuningHz = hz;
        }
    };
}
