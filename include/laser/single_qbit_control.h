#pragma once
#include <functional>

#include "quantum_computer/time_master.h"
#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "pulse_command.h"
#include "two_photon_laser.h"

#include "quantum_gates/gate_diag.h"

namespace KetCat
{
    /// @brief Single-qubit quantum control layer for Raman/STIRAP-driven operations.
    ///
    /// @details
    ///   This class translates abstract logical pulse commands into
    ///   physically realizable laser-driven dynamics and propagates
    ///   the quantum state under the resulting time-dependent Hamiltonian.
    ///
    ///   The controller operates on a fixed neutral-atom manifold and
    ///   supports single-qubit rotations implemented through:
    ///
    ///     • Raman-driven X/Y rotations
    ///     • Virtual Z rotations via rotating-frame updates
    ///     • STIRAP / Inverted STIRAP pulse sequencing
    ///
    ///   Main responsibilities:
    ///
    ///     • Interpret logical PulseCommand operations
    ///     • Configure two-photon laser interactions
    ///     • Construct and update the driven Hamiltonian
    ///     • Propagate the wavefunction using Crank–Nicolson evolution
    ///     • Maintain coherent rotating-frame phase tracking
    ///
    ///   Physical model:
    ///
    ///     • Logical qubit states are embedded in a three-level ladder system
    ///     • X/Y rotations are realized through coherent Raman coupling
    ///     • Z rotations are implemented virtually without population transfer
    ///     • Laser phases encode the effective rotation axis
    ///
    ///   Numerical evolution:
    ///
    ///     H(t) Ψ(t) = i ∂Ψ(t)/∂t
    ///
    ///   where:
    ///
    ///     • H(t) is the time-dependent driven Hamiltonian
    ///     • Ψ(t) is the propagated quantum state
    ///
    ///   The propagation is performed incrementally using a Crank–Nicolson solver.
    ///
    ///   STIRAP sequencing:
    ///
    ///     • Successive rotations alternate between:
    ///         - STIRAP
    ///         - Inverted STIRAP
    ///
    ///   based on the accumulated total rotation angle:
    ///
    ///     θ_accum ∈ [0, 2π)
    ///
    /// @tparam Config
    ///   Neutral atom configuration describing:
    ///     • energy levels
    ///     • dipole couplings
    ///     • logical qubit mapping
    template <NeutralAtomTypeConfig Config>
    class SingleQubitControl
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;
        using OperationHilbertSpace = typename Manifold::SingleAtomOperationHilbertSpace;

        /// @brief Atomic energy levels in Hartree atomic units.
        ///
        /// @details
        ///   Retrieved from the configured neutral-atom manifold.
        const std::array<real_t, ConfigType::LevelCount> m_energies =
            Manifold::getHartreeEnergies();

        /// @brief Electric dipole coupling matrix between atomic states.
        ///
        /// @details
        ///   Matrix elements define laser-induced transition strengths:
        ///
        ///     μᵢⱼ = ⟨i| d̂ |j⟩
        const matrix_t<ConfigType::LevelCount> m_dipoleMatrix =
            Manifold::getDipoleMatrix();

        /// @brief Peak Rabi frequency used for generated pulses (Hz).
        ///
        /// @details
        ///   Defines the maximum laser coupling amplitude:
        ///
        ///     Ω_max / (2π)
        real_t m_peakRabiHz = 100e6;

        /// @brief Common Raman detuning from the intermediate state (Hz).
        ///
        /// @details
        ///   Applied symmetrically to the two-photon interaction.
        real_t m_commonDetuningHz = 0;

        /// Tracks the effective logical qubit frame phase.
        ///
        /// Includes:
        ///   • Virtual Z gate updates
        ///   • Pulse-induced residual phase shifts
        ///   • Calibration frame corrections
        ///
        /// Subsequent Raman pulse phases are expressed relative
        /// to this accumulated logical frame.
        real_t m_framePhase = 0.0;

        /// @brief Accumulated STIRAP sequencing angle.
        ///
        /// @details
        ///   Used to alternate between:
        ///
        ///     • STIRAP
        ///     • Inverted STIRAP
        ///
        ///   based on the accumulated total rotation angle:
        ///
        ///     θ_accum ∈ [0, 2π)
        static inline real_t m_stirapCycleTheta = 0.0;

        /// @brief Select STIRAP protocol based on accumulated rotation angle.
        ///
        /// @param targetTheta
        ///   Requested logical rotation angle in radians.
        ///
        /// @return
        ///   Selected two-photon transfer protocol.
        ///
        /// @details
        ///   Protocol selection alternates between:
        ///
        ///     • STIRAP
        ///     • Inverted STIRAP
        ///
        ///   over successive operations in order to balance
        ///   adiabatic transfer directionality.
        ///
        ///   The accumulated angle is wrapped periodically:
        ///
        ///     θ_accum → θ_accum mod 2π
        static TwoPhotonProtocol handleStirapTheta(real_t targetTheta)
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

            /// Wrap accumulated angle into [0, 2π)
            if (ConstexprMath::floatNear(
                m_stirapCycleTheta,
                2 * ConstexprMath::Pi))
            {
                m_stirapCycleTheta = 0.0;
            }

            return Result;
        }

        /// @brief Construct physical laser configuration for a pulse command.
        ///
        /// @param command
        ///   Logical pulse command to translate.
        ///
        /// @return
        ///   Fully configured two-photon laser interaction.
        ///
        /// @details
        ///   Rotation-axis encoding:
        ///
        ///     • X rotation → phase = 0
        ///     • Y rotation → phase = π/2
        ///
        ///   The final physical laser phase additionally includes
        ///   the accumulated rotating-frame correction:
        ///
        ///     φ_total = φ_axis - φ_frame
        ///
        ///   The generated configuration includes:
        ///
        ///     • Energy levels
        ///     • Dipole couplings
        ///     • Rabi frequencies
        ///     • Detuning
        ///     • Pulse phases
        ///     • STIRAP protocol selection
        TwoPhotonConfig prepareLaserConfig(const PulseCommand& command)
        {
            /// Determine phase corresponding to the requested rotation axis
            real_t AxisPhase =
                (command.m_axis == RotationAxis::Y)
                ? (ConstexprMath::Pi / 2.0)
                : 0.0;

            /// Include accumulated rotating-frame correction - TODO cleanup and measure actual gate phase shifts
            real_t TotalLaserPhase = AxisPhase - m_framePhase;
            
            m_framePhase += AxisPhase;

            /// Configure two-photon interaction
            TwoPhotonConfig LaserConfig;

            LaserConfig.m_Level1Energy =
                m_energies[ConfigType::Logical0Level];

            LaserConfig.m_Level2Energy =
                m_energies[ConfigType::Logical0Level + 1];

            LaserConfig.m_Level3Energy =
                m_energies[ConfigType::Logical1Level];

            LaserConfig.m_Mu12 =
                m_dipoleMatrix
                [ConfigType::Logical0Level]
                [ConfigType::Logical0Level + 1]
                .re;

            LaserConfig.m_Mu23 =
                m_dipoleMatrix
                [ConfigType::Logical0Level + 1]
                [ConfigType::Logical1Level]
                .re;

            LaserConfig.m_peakRabiFrequency =
                Units::omegaAuFromHz(m_peakRabiHz);

            LaserConfig.m_commonDetuning =
                Units::omegaAuFromHz(m_commonDetuningHz);

            LaserConfig.m_targetTheta =
                command.m_rotationAngleRad;

            LaserConfig.m_pumpPhase =
                TotalLaserPhase;

            LaserConfig.m_protocol =
                handleStirapTheta(command.m_rotationAngleRad);

            std::cout
                << "Laser parameters:"
                << "\n  Peak Rabi Frequency (Hz): " << m_peakRabiHz
                << "\n  Common Detuning (Hz): " << m_commonDetuningHz
                << "\n  Target Theta (rad): " << command.m_rotationAngleRad
                << "\n  Total Laser Phase (rad): " << m_framePhase
                << "\n  Protocol: "
                << (LaserConfig.m_protocol ==
                    TwoPhotonProtocol::STIRAP
                    ? "STIRAP"
                    : "Inverted STIRAP")
                << std::endl;

            return LaserConfig;
        }

    public:
        using ControlLaserArray = std::array<LaserPulse, ConfigType::LevelCount - 1>;

        /// @brief Callback type for observing pulse evolution.
        using CallbackType = std::function<void(const ControlLaserArray&)>;

        /// @brief Apply a logical pulse command to the quantum state.
        ///
        /// @param command
        ///   Requested logical rotation.
        ///
        /// @param psi
        ///   Input/output quantum state vector.
        ///
        /// @param callback
        ///   Optional evolution callback invoked after each propagation step.
        ///
        /// @details
        ///   Supported operations:
        ///
        ///     • Z rotations:
        ///
        ///         Implemented virtually by updating the rotating frame:
        ///
        ///             φ_frame ← φ_frame + θ
        ///
        ///         No physical laser pulse or propagation is performed.
        ///
        ///     • X/Y rotations:
        ///
        ///         1. Generate Raman/STIRAP pulse envelopes
        ///         2. Construct driven Hamiltonian
        ///         3. Propagate the wavefunction incrementally
        ///         4. Emit intermediate callback updates
        ///
        ///   Time evolution:
        ///
        ///     Ψ(t + Δt) = U(Δt) Ψ(t)
        ///
        ///   where U(Δt) is computed using the Crank–Nicolson method.
        void applyPulseCommand(
            const PulseCommand& command,
            StateVector<OperationHilbertSpace>& psi,
            const CallbackType& callback)
        {
            if (command.m_axis == RotationAxis::Z)
            {
                /// Update rotating-frame phase
                m_framePhase += command.m_rotationAngleRad;

                /// Virtual Z requires no physical pulse
                return;
            }

            TwoPhotonLaser LaserEnvelope(
                prepareLaserConfig(command));

            TimeMaster::Clock().setInstructionTimeLimit(
                LaserEnvelope.getFullTransferTime());

            real_t GateTime =
                LaserEnvelope.getTransitionTimeLimit();
            
            while (TimeMaster::Clock().getCurrentInstructionTime() < GateTime)
            {
                auto [Pump, Stokes] =
                    LaserEnvelope(
                        TimeMaster::Clock()
                        .getCurrentStirapCycleTime());

                ControlLaserArray Lasers;
                Lasers[ConfigType::Logical0Level] = Pump;
                Lasers[ConfigType::Logical0Level + 1] = Stokes;

                callback(Lasers);

                TimeMaster::Clock().tick();
            }

            /// Reset instruction-local timing state
            TimeMaster::Clock().resetCurrentInstructionClock();

        }

        /// @brief Set peak Rabi frequency for generated pulses.
        ///
        /// @param hz
        ///   Peak Rabi frequency in Hz.
        ///
        /// @details
        ///   Controls the maximum Raman coupling strength:
        ///
        ///     Ω_max / (2π)
        void setPeakRabiHz(real_t hz)
        {
            m_peakRabiHz = hz;
        }

        /// @brief Set common Raman detuning.
        ///
        /// @param hz
        ///   Detuning in Hz.
        ///
        /// @details
        ///   Applied uniformly to the intermediate-state detuning:
        ///
        ///     Δ = ω_laser - ω_transition
        void setDetuningHz(real_t hz)
        {
            m_commonDetuningHz = hz;
        }
    };
}
