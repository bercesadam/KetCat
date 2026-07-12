#pragma once
#include <optional>

#include "atomic_units.h"
#include "two_photon_laser.h"
#include "quantum_processor/time_master.h"
#include "compiler/physical_instruction.h"
#include "local_space/neutral_atom_manifold.h"

#include "codegen/gate_calibration_map.h"


namespace KetCat
{
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

        /// @brief Hartree energies for the relevant atomic levels.
        const std::array<real_t, ConfigType::LevelCount> m_energies =
            Manifold::getHartreeEnergies();
        
		/// @brief Electric dipole transition matrix between eigenstates.
        const square_matrix_t<ConfigType::LevelCount> m_dipoleMatrix =
            Manifold::getDipoleMatrix();

		/// @brief Peak Rabi frequency Ω₀ / (2π) for the Raman transition.
        real_t m_peakRabiHz{};

		/// @brief Common detuning Δ / (2π) for the intermediate state in the Raman transition.
        real_t m_commonDetuningHz{};

		/// @brief Per-atom phases for the rotating frame, used to track virtual Z rotations and maintain coherence.
        std::array<real_t, AtomCount> m_RWAFramePhases{};


        /// @brief Construct physical laser configuration for a pulse command.
        ///
        /// @details
        ///    Calculates the required laser frequencies and phases.
        ///    ωp = (E₁ - E₀) + Δ
        ///    ωs = (E₂ - E₁) - Δ
        ///
        ///    The rotating frame update ensures phase coherence:
        ///      φ_pump = φ_instr - φ_frame
        TwoPhotonConfig prepareLaserConfig(const PhysicalInstruction& instruction)
        {
            TwoPhotonConfig LaserConfig;

            LaserConfig.m_GroundLevelIndex = (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
				? ConfigType::Logical1Level
				: ConfigType::Logical0Level;

            real_t PeakRabiHzP = 0.0;
            real_t PeakRabiHzS = 0.0;
            real_t CommonDetuningHz = 0.0;

            if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
            {
                PeakRabiHzP = 500e6;
                PeakRabiHzS = PeakRabiHzP;
                CommonDetuningHz = 10000e6;
            }
            else
            {
                PeakRabiHzP = 150e6;
                PeakRabiHzS = PeakRabiHzP;
                CommonDetuningHz = 1500e6;
            }

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

            LaserConfig.m_peakRabiFrequencyP = Units::omegaAuFromHz(PeakRabiHzP);
            LaserConfig.m_peakRabiFrequencyS = Units::omegaAuFromHz(PeakRabiHzS);
            LaserConfig.m_commonDetuning = Units::omegaAuFromHz(CommonDetuningHz);

            LaserConfig.m_pumpPhase = instruction.m_phases[0].m_phase;
            LaserConfig.m_pumpPhase2 = instruction.m_phases[1].m_phase;
            LaserConfig.m_pumpPhase2Timing = instruction.m_phases[1].m_time;
            LaserConfig.m_protocol = TwoPhotonProtocol::Simultaneous;

            LaserConfig.m_targetTheta = instruction.m_theta;
          
            return LaserConfig;
        }

    public:
        constexpr void initializeAtomAsLogical1(natural_t index)
        {
			(void)index; // Placeholder for any future logic related to initializing the frame phase for logical |1⟩ states.
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
        std::optional<TwoPhotonLaserEnvelope> calculateLaserEnvelope(PhysicalInstruction& instruction)
        {
            const natural_t targetIndex = instruction.m_targets[0];
			real_t& framePhase = m_RWAFramePhases[targetIndex];

            if (targetIndex >= AtomCount)
            {
                return std::nullopt;
            }

			std::cout << "Rotating frame phase before instruction: " << framePhase << " radians" << std::endl;

            // Handle Bloch-Z rotations (Virtual Z)
            if (instruction.m_type == PhysicalInstructionType::VirtualZ)
            {
                framePhase += instruction.m_theta;
                return std::nullopt;
            }
           
            // Apply rotating frame correction: φ_eff = φ_laser - φ_frame
            for (TimeDependentPhase& p : instruction.m_phases)
            {
                p.m_phase -= framePhase;
            }

            // Generate physical laser parameters via Builder
            TwoPhotonPulseBuilder LaserBuilder(prepareLaserConfig(instruction));
            TwoPhotonLaserEnvelope LaserEnvelope = LaserBuilder.build();

            // Apply post-gate phase error RWA frame corrections, in case the gates are calibrated
            applyPostGatePhaseCorrection(instruction);

            return LaserEnvelope;
        }

    private:
        /// @brief Applies dynamic phase corrections to the Rotating Wave Approximation (RWA) frames following entangling or local operations.
        ///
        /// @details
        /// When `CALIBRATED_GATES` is active and `DISABLE_PHASE_GATE_CORRECTION` is omitted, this function injects 
        /// compensatory phase updates directly into the RWA frame tracking registers of the target atoms. 
        /// For Rydberg excitations (e.g., Controlled-Phase gates), fixed control and target frame errors are added. 
        /// For Raman rotations, a nonlinear phase error is evaluated dynamically via interpolation as a function of the rotation angle.
        ///
        /// @param instruction The executed physical control pulse structure containing target indices and operational parameters.
        constexpr void applyPostGatePhaseCorrection(PhysicalInstruction& instruction) noexcept
        {
#if defined(CALIBRATED_GATES) && !defined(DISABLE_PHASE_GATE_CORRECTION) 

                if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
                {
                    static constexpr TwoQubitCalibResult CzError = GateCalibrationTable::getCPhaseCalib();

					// 1. Apply the local Stark-shift phase errors to the control and target atoms
                    constexpr real_t deltaCzPhase = CzError.m_actualCzPhase - ConstexprMath::Pi;

					// 2. Distribute the collective phase error symmetrically between the two atoms
                    constexpr real_t collectiveCorrection = deltaCzPhase / 2.0;

					// Control atom receives its local Stark-shift error AND half of the collective error
                    m_RWAFramePhases[instruction.m_targets[0]] += (CzError.m_controlFramePhaseError + collectiveCorrection);

					// Target atom receives its local Stark-shift error AND half of the collective error
                    m_RWAFramePhases[instruction.m_targets[1]] += (CzError.m_targetFramePhaseError + collectiveCorrection);

                    std::cout << "CPHASE RWA frame correction:" << std::endl;
                    std::cout << "  Measured Pure CZ Phase: " << CzError.m_actualCzPhase << " rad (Error: " << deltaCzPhase << " rad)" << std::endl;
                    std::cout << "  Applied Symmetric Collective Correction: +" << collectiveCorrection << " rad/atom" << std::endl;
                    std::cout << "  Total Control Frame Shift: " << (CzError.m_controlFramePhaseError + collectiveCorrection) << " rad" << std::endl;
                    std::cout << "  Total Target Frame Shift: " << (CzError.m_targetFramePhaseError + collectiveCorrection) << " rad" << std::endl;
                }
            else if (instruction.m_type == PhysicalInstructionType::RamanRotation)
            {
                static auto PhaseInterpolator = GateCalibrationTable::getRxCalib();
                real_t PhaseError = PhaseInterpolator.evaluate(instruction.m_theta);
                m_RWAFramePhases[instruction.m_targets[0]] += PhaseError;
                std::cout << "Applied RWA frame correction after RamanRotation for theta " << instruction.m_theta << ": +" << PhaseError << " rad" << std::endl;
            }
#endif
        }
    };
}
