#pragma once
#include <functional>

#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "pulse_command.h"
#include "stirap.h"

namespace KetCat
{
	constexpr real_t TimeStep = 50; //AU

	/// @brief Protocol for single qubit operations using a three-level STIRAP scheme.
	template <NeutralAtomTypeConfig Config>
	class SingleQubitProtocol
	{
		using ConfigType = std::remove_cvref_t<decltype(Config)>;
		using Manifold = NeutralAtomManifold<Config>;
		using OperationHilbertSpace = Manifold::SingleAtomOperationHilbertSpace;

		const std::array<real_t, ConfigType::LevelCount> Energies = Manifold::getHartreeEnergies();
		const matrix_t<ConfigType::LevelCount> DipoleMatrix = Manifold::getDipoleMatrix();

		using CallbackType = std::function<void(real_t, const StateVector<OperationHilbertSpace>&, const LaserPulse&, const LaserPulse&)>;

 public:
		void applyPulseCommand(const PulseCommand& command, StateVector<OperationHilbertSpace>& psi, const CallbackType& callback)
		{
			if (command.m_axis == RotationAxis::X || command.m_axis == RotationAxis::Y)
			{
				STIRAPConfig LaserConfig;
				LaserConfig.m_Level1Energy = Energies[ConfigType::Logical0Level];
				LaserConfig.m_Level2Energy = Energies[ConfigType::Logical0Level + 1];
				LaserConfig.m_Level3Energy = Energies[ConfigType::Logical1Level];
				LaserConfig.m_Mu12 = DipoleMatrix[ConfigType::Logical0Level][ConfigType::Logical0Level + 1].re;
				LaserConfig.m_Mu23 = DipoleMatrix[ConfigType::Logical0Level + 1][ConfigType::Logical1Level].re;
				LaserConfig.m_pumpPhase = (command.m_axis == RotationAxis::X) ? 0.0 : ConstexprMath::Pi / 2.0; // Phase for X or Y rotation.
				LaserConfig.m_rabiFrequency_Hz = 1E9; // Target Rabi frequency in Hz (can be adjusted based on experimental parameters).
				LaserConfig.m_stokesDetuning = 0.0; // On-resonance for the Stokes pulse.
				LaserConfig.m_targetTheta = command.m_rotationAngleRad; // Target rotation angle in radians.

				STIRAPLaser<Config> Laser(LaserConfig);

				real_t GateTime = Laser.getTransitionTimeLimit();
				real_t CurrentTime = 0.0;

				LaserPulse Pump, Stokes;

				while (CurrentTime < GateTime)
				{
					std::tie(Pump, Stokes) = Laser(CurrentTime);
					std::array<LaserPulse, ConfigType::LevelCount - 1> Lasers = { Pump, Stokes };
					MultiRwaRabiHamiltonian<ConfigType::LevelCount> H(Energies, DipoleMatrix, Lasers);

					CrankNicolsonSolver<OperationHilbertSpace> solver(H.getMatrix(), TimeStepAu);
					psi = solver(psi);

					callback(CurrentTime, psi, Pump, Stokes);

					CurrentTime += TimeStep;
				}
				callback(CurrentTime, psi, Pump, Stokes);
			}
		}
	};	
}
