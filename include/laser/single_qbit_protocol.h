#pragma once
#include <functional>

#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "pulse_command.h"
#include "stirap.h"

namespace KetCat
{
    /**
     * @brief Translates logical pulse commands into physical time-evolution.
     * @details Implements physical X/Y rotations via TDSE and Virtual Z gates via phase shifts.
     */
    template <NeutralAtomTypeConfig Config>
    class SingleQubitProtocol
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;
        using OperationHilbertSpace = typename Manifold::SingleAtomOperationHilbertSpace;

        const std::array<real_t, ConfigType::LevelCount> m_energies = Manifold::getHartreeEnergies();
        const matrix_t<ConfigType::LevelCount> m_dipoleMatrix = Manifold::getDipoleMatrix();

        real_t m_timeStepAu = 0.1;
        real_t m_peakRabiHz = 50e6;
        real_t m_commonDetuningHz = 2e9;

        // This tracks the accumulated phase for Virtual Z gates to be applied 
        // to subsequent physical X/Y pulses.
        real_t m_framePhase = 0.0;

    public:
        using CallbackType = std::function<void(real_t, const StateVector<OperationHilbertSpace>&, const LaserPulse&, const LaserPulse&, const bool)>;

        SingleQubitProtocol(real_t timeStepAu = 0.1) : m_timeStepAu(timeStepAu) {}

        void applyPulseCommand(const PulseCommand& command, StateVector<OperationHilbertSpace>& psi, const CallbackType& callback)
        {
            if (command.m_axis == RotationAxis::Z)
            {
				// 1. Update the accumulated frame phase for virtual Z rotation
                m_framePhase += command.m_rotationAngleRad;

				// 2. Apply the phase shift to the logical |1⟩ level in the state vector.
                psi[ConfigType::Logical1Level] = psi[ConfigType::Logical1Level] * complex_t::fromPolar(1.0, command.m_rotationAngleRad);

                // No laser pulse needed
                return;
            }

            // --- PHYSICAL X / Y GATES ---

			// For X/Y rotations, we generate a two-photon Raman pulse sequence. The phase of the drive is determined by the target rotation axis:
            real_t axisPhase = (command.m_axis == RotationAxis::Y) ? (ConstexprMath::Pi / 2.0) : 0.0;
            real_t totalLaserPhase = axisPhase + m_framePhase;

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
            laserConfig.m_protocol = TwoPhotonPulseProtocol::Simultaneous;

            GenericTwoPhotonLaser laser(laserConfig);

            real_t gateTime = laser.getTransitionTimeLimit();
            real_t currentTime = 0.0;

			auto [pump, stokes] = laser(0);
            std::array<LaserPulse, ConfigType::LevelCount - 1> lasers;
            lasers[ConfigType::Logical0Level] = pump;
            lasers[ConfigType::Logical0Level + 1] = stokes;

            MultiRwaRabiHamiltonian<ConfigType::LevelCount> H(m_energies, m_dipoleMatrix, lasers);

            while (currentTime < gateTime)
            {
                std::tie(pump, stokes) = laser(currentTime);
                
                lasers[ConfigType::Logical0Level] = pump;
                lasers[ConfigType::Logical0Level + 1] = stokes;

                H.updateOffDiagonal(lasers);
                CrankNicolsonSolver<OperationHilbertSpace> solver(H.getMatrix(), m_timeStepAu);
                
                psi = solver(psi);
                callback(currentTime, psi, pump, stokes, false);

                currentTime += m_timeStepAu;
            }

            callback(currentTime, psi, pump, stokes, true);
        }

        void setPeakRabiHz(real_t hz) { m_peakRabiHz = hz; }
        void setDetuningHz(real_t hz) { m_commonDetuningHz = hz; }
    };
}
