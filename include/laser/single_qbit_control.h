#pragma once
#include <functional>

#include "systems/time_master.h"
#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "pulse_command.h"
#include "two_photon_laser.h"


namespace KetCat
{
    /// @brief Control layer that maps logical pulse commands to physical time evolution.
    ///
    /// @details
    ///   This class translates high-level single-qubit rotation commands into
    ///   laser-driven dynamics and propagates the quantum state accordingly.
    ///
    ///   Responsibilities:
    ///
    ///     • Interpret PulseCommand (X, Y, Z rotations)
    ///     • Convert commands into two-photon Raman / STIRAP laser parameters
    ///     • Build and update the time-dependent Hamiltonian
    ///     • Execute time evolution using a numerical solver (Crank–Nicolson)
    ///     • Maintain a rotating frame for Virtual Z gates
    ///
    ///   Gate implementation:
    ///
    ///     • X / Y gates:
    ///         Implemented as physical Raman transitions using dynamically
    ///         generated laser pulses. The rotation axis is encoded as a phase
    ///         of the driving field.
    ///
    ///     • Z gates (virtual):
    ///         Implemented as a frame update (phase accumulation) without
    ///         applying any physical pulse. The phase is applied directly to
    ///         the state and tracked for subsequent operations.
    ///
    ///   Execution model:
    ///
    ///     • A pulse command is expanded into a time-dependent laser sequence
    ///     • The Hamiltonian is updated at each time step
    ///     • The state is propagated incrementally
    ///     • Optional callback provides access to intermediate states and pulses
    ///
    ///   Notes:
    ///
    ///     • This is a single-qubit control layer operating on a fixed atomic manifold
    ///     • The implementation assumes a three-level system for Raman transitions
    ///     • STIRAP / Inverted STIRAP protocols are alternated between calls
    ///
    /// @tparam Config
    ///   Neutral atom configuration describing level structure and couplings.
    template <NeutralAtomTypeConfig Config>
    class SingleQubitControl
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;
        using OperationHilbertSpace = typename Manifold::SingleAtomOperationHilbertSpace;

        /// @brief Energy levels (Hartree units) for the configured atom.
        const std::array<real_t, ConfigType::LevelCount> m_energies = Manifold::getHartreeEnergies();

        /// @brief Dipole coupling matrix between atomic levels.
        const matrix_t<ConfigType::LevelCount> m_dipoleMatrix = Manifold::getDipoleMatrix();

        /// @brief Peak Rabi frequency used for generated pulses (Hz).
        real_t m_peakRabiHz = 100e6;

        /// @brief Common detuning applied to Raman transitions (Hz).
        real_t m_commonDetuningHz = 0;

        /// @brief Accumulated phase of the rotating frame for Virtual Z gates.
        ///
        /// @details
        ///   This phase is added to subsequent physical pulses to ensure
        ///   correct phase coherence across operations.
        real_t m_framePhase = 0.0;

    public:
        /// @brief Callback signature for monitoring time evolution.
        ///
        /// @param time            Current simulation time
        /// @param state           Current quantum state
        /// @param pump            Pump laser pulse
        /// @param stokes          Stokes laser pulse
        /// @param isFinalStep     True if this is the last step of the pulse
        using CallbackType = std::function<void(const StateVector<OperationHilbertSpace>&, const LaserPulse&, const LaserPulse&, const bool)>;

        /// @brief Apply a single pulse command to the quantum state.
        ///
        /// @param command
        ///   Logical rotation command (axis + angle).
        ///
        /// @param psi
        ///   Input/output state vector. Updated in-place during evolution.
        ///
        /// @param callback
        ///   User callback invoked at each time step:
        ///     (time, state, pump pulse, stokes pulse, isFinalStep)
        ///
        /// @details
        ///   Behavior depends on the rotation axis:
        ///
        ///     • Z:
        ///         Applies a virtual phase shift to the logical |1⟩ state and updates
        ///         the internal frame phase. No time evolution is performed.
        ///
        ///     • X / Y:
        ///         Generates a Raman pulse sequence, constructs the Hamiltonian,
        ///         and propagates the state over the pulse duration.
        void applyPulseCommand(const PulseCommand& command, StateVector<OperationHilbertSpace>& psi, const CallbackType& callback)
        {
            if (command.m_axis == RotationAxis::Z)
            {
                /// Update accumulated frame phase for virtual Z rotation
                m_framePhase += command.m_rotationAngleRad;

                /// No physical pulse required
                return;
            }

            /// Determine phase corresponding to rotation axis
            real_t axisPhase = (command.m_axis == RotationAxis::Y) ? (ConstexprMath::Pi / 2.0) : 0.0;

			// Total laser phase includes both the desired rotation axis and the accumulated frame phase
            real_t totalLaserPhase = axisPhase - m_framePhase;

            /// Configure two-photon Raman interaction
            TwoPhotonConfig laserConfig;
            laserConfig.m_Level1Energy = m_energies[ConfigType::Logical0Level];
            laserConfig.m_Level2Energy = m_energies[ConfigType::Logical0Level + 1];
            laserConfig.m_Level3Energy = m_energies[ConfigType::Logical1Level];

            laserConfig.m_Mu12 = m_dipoleMatrix[ConfigType::Logical0Level][ConfigType::Logical0Level + 1].re;
            laserConfig.m_Mu23 = m_dipoleMatrix[ConfigType::Logical0Level + 1][ConfigType::Logical1Level].re;

            laserConfig.m_peakRabiFrequency = Units::omegaAuFromHz(m_peakRabiHz);
            laserConfig.m_commonDetuning = Units::omegaAuFromHz(m_commonDetuningHz);
            laserConfig.m_targetTheta = command.m_rotationAngleRad;
            laserConfig.m_pumpPhase = totalLaserPhase;

			static bool wasFractionalStirap = false;

            /// Alternate between STIRAP and inverted STIRAP
            static bool even = false;
            laserConfig.m_protocol = even ? TwoPhotonProtocol::InvertedSTIRAP : TwoPhotonProtocol::STIRAP;

            if (command.m_rotationAngleRad >= ConstexprMath::Pi - 1e-7)
            {
                wasFractionalStirap = false;
            }
            else
            {   
                wasFractionalStirap = true;
            }

            if (!wasFractionalStirap)
            {
                even = !even;
            }


            std::cout << "Laser parameters:" << "\n  Peak Rabi Frequency (Hz): " << m_peakRabiHz
                << "\n  Common Detuning (Hz): " << m_commonDetuningHz
                << "\n  Target Theta (rad): " << command.m_rotationAngleRad
                << "\n  Total Laser Phase (rad): " << totalLaserPhase
                << "\n  Protocol: " << (laserConfig.m_protocol == TwoPhotonProtocol::STIRAP ? "STIRAP" : "Inverted STIRAP")
				<< "\n  Fractional STIRAP: " << (wasFractionalStirap ? "Yes" : "No")
                << std::endl;

            TwoPhotonLaser laser(laserConfig);

            real_t gateTime = laser.getTransitionTimeLimit();

			std::cout << "Gate time (s): " << gateTime << std::endl;
			std::cout << "Current instruction time: " << TimeMaster::Clock().getCurrentInstructionTime() << std::endl;

            auto [pump, stokes] = laser(0);
            std::array<LaserPulse, ConfigType::LevelCount - 1> lasers;
            lasers[ConfigType::Logical0Level] = pump;
            lasers[ConfigType::Logical0Level + 1] = stokes;

            MultiRwaRabiHamiltonian<ConfigType::LevelCount> Hamiltonian(m_energies, m_dipoleMatrix, lasers);

            while (TimeMaster::Clock().getCurrentInstructionTime() < gateTime)
            {
                std::tie(pump, stokes) = laser(TimeMaster::Clock().getCurrentInstructionTime());

                lasers[ConfigType::Logical0Level] = pump;
                lasers[ConfigType::Logical0Level + 1] = stokes;

                Hamiltonian.updateOffDiagonal(lasers);
                Hamiltonian.updateMainDiagonal(lasers);

                CrankNicolsonSolver<OperationHilbertSpace> solver(Hamiltonian.getMatrix(), TimeMaster::Clock().getTimeStep());

                psi = solver(psi);
                callback(psi, pump, stokes, false);

                TimeMaster::Clock().tick();
            }

			// Ensure final callback is invoked with the last state and pulse parameters, marking the end of the pulse sequence
            callback(psi, pump, stokes, true);

            if (command.m_rotationAngleRad >= ConstexprMath::Pi - 1e-7)
            {
                TimeMaster::Clock().resetCurrentInstructionClock();
            }
            else
            {
				std::cout << "valami" << std::endl;    
            }
        }

        /// @brief Set peak Rabi frequency used for generated pulses (Hz).
        void setPeakRabiHz(real_t hz)
        { 
            m_peakRabiHz = hz;
        }

        /// @brief Set common detuning for Raman transitions (Hz).
        void setDetuningHz(real_t hz)
        {
            m_commonDetuningHz = hz;
        }
    };
}