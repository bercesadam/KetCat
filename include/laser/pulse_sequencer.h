#pragma once
#include <optional>
#include "two_photon_laser.h"
#include "quantum_processor/time_master.h"
#include "compiler/physical_instruction.h"
#include "operation_space/neutral_atom_manifold.h"

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

        const std::array<real_t, ConfigType::LevelCount> m_energies =
            Manifold::getHartreeEnergies();

        const square_matrix_t<ConfigType::LevelCount> m_dipoleMatrix =
            Manifold::getDipoleMatrix();

        real_t m_peakRabiHz = 50e6;
        real_t m_commonDetuningHz = 500e6;

        std::array<real_t, AtomCount> m_RWAFramePhases;


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

            if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
            {
                m_peakRabiHz = 500e6;
                m_commonDetuningHz = 1000e6;
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

            LaserConfig.m_peakRabiFrequency = Units::omegaAuFromHz(m_peakRabiHz);
            LaserConfig.m_commonDetuning = Units::omegaAuFromHz(m_commonDetuningHz);

            LaserConfig.m_pumpPhase = instruction.m_phase;
            LaserConfig.m_protocol = TwoPhotonProtocol::Simultaneous;
          
            return LaserConfig;
        }

    public:
        constexpr void initializeAtomAsLogical1(natural_t index)
        {
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
			real_t& framePhase = m_RWAFramePhases[targetIndex];

            if (targetIndex >= AtomCount)
            {
                return std::nullopt;
            }

            // Handle Bloch-Z rotations (Virtual Z)
            if (instruction.m_type == PhysicalInstructionType::VirtualZ)
            {
                framePhase += instruction.m_theta;
                return std::nullopt;
            }
           
            // Apply rotating frame correction: φ_eff = φ_laser - φ_frame
            instruction.m_phase -= framePhase;

            // Generate physical laser parameters via Builder
            TwoPhotonPulseBuilder LaserBuilder(prepareLaserConfig(instruction));
            TwoPhotonLaserEnvelope LaserEnvelope = LaserBuilder.build();

            // For X/Y rotations, the instruction phase defines the axis
            // Update frame phase if relative phase tracking is required
            framePhase += instruction.m_phase;

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
