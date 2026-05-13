
#pragma once

#include <array>
#include <variant>
#include <tuple>
#include <utility>
#include <type_traits>

#include "operation_space/neutral_atom_manifold.h"
#include "operation_space/utils/subspace_operations.h"

#include "quantum_circuit.h"
#include "compiler/gate_compiler.h"

#include "laser/pulse_sequencer.h"

#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"


namespace KetCat
{
    //============================================================
    // Observer system
    //============================================================

    template<typename StateType>
    class ExecutionObserver
    {
    public:

        virtual ~ExecutionObserver() = default;

        virtual void onTimeStep(
            real_t globalTime,
            const StateType& psi) = 0;
    };


    //============================================================
    // QPU
    //============================================================

    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class QuantumProcessor
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;

        StateVector<typename GlobalStateManager::FullHilbertSpace> m_GlobalStateVector;

        NeutralAtomManifold<Config> m_Manifold;

        LaserPulseSequencer<Config, QubitCount> m_laserSequencer;

    public:

        constexpr QuantumProcessor(real_t dt)
        {
            TimeMaster::Clock().init(dt);
            StateVector<decltype(m_Manifold)::SingleAtomOperationHilbertSpace> OneAtomSeed =
                m_Manifold.getOperationSeed();
            m_GlobalStateVector = GlobalStateManager::productStateFromSeed(OneAtomSeed);
        }

        // Iterate thru the logical gates and execute their associated physical actions
        template<natural_t OpCount>
        void execute(const QuantumCircuit<QubitCount, OpCount>& circuit)
        {
           
            for (const auto& op : circuit.operations())
            {
                std::visit(
                    [&](const auto& gate)
                {
                    executeGate(gate);
                }, op);
            }
        }

    private:
        // Translate the logical gate into a series of physical actions
        template<typename GateOp>
        void executeGate(const GateOp& gate)
        {
            
            auto [ PhysicalInstructions, InstructionCount]
                = GateCompiler::compile(gate);

            for (natural_t i = 0; i < InstructionCount; ++i)
            {
                executeInstruction(PhysicalInstructions[i]);
            }
        }

        // Translate the physical instruction into a series of laser pulses
        void executeInstruction(const PhysicalInstruction& instruction)
        {
            
            auto PulseEnvelope = m_laserSequencer.calculateLaserEnvelope(instruction);

            // In case of Virtual Z pulses, only internal changes are taken
            // in the laser controller and we don't need TDSE evolution
            if (!PulseEnvelope)
            {
                return;
            }

            TwoPhotonLaserEnvelope& Envelope = *PulseEnvelope;
            real_t TransitionTimeLimit = Envelope.getTransitionTimeLimit();
            real_t TimeShift = Envelope.setStartTime(); // if adiabatic continuity needs to be ensured

            while (TimeMaster::Clock().getCurrentInstructionTime() < TransitionTimeLimit)
            {
                // TODO probably for handling Rydberg excitations we need a better laser pulse structure
                // not getting along by Pump and Stokes for the lower 3 states
                auto [Pump, Stokes] =
                    Envelope(TimeMaster::Clock().getCurrentInstructionTime() + TimeShift);

                // Compile array with the drive laser for each energy level of the operational space
                typename MultiRwaRabiHamiltonian<ConfigType::LevelCount>::template laser_array_t Lasers;
                Lasers[ConfigType::Logical0Level] = Pump;
                Lasers[ConfigType::Logical0Level + 1] = Stokes;

                // Perform one time step of the Schrödinger evolution
                evolveGlobalState(Lasers, instruction.m_targets[0]);

                // observer

                TimeMaster::Clock().tick();
            }

            /// Reset instruction-local timing state
            TimeMaster::Clock().resetCurrentInstructionClock();
        }

        void evolveGlobalState(const MultiRwaRabiHamiltonian<ConfigType::LevelCount>::laser_array_t lasers, const natural_t affectedQubit)
        {
            static const std::array<real_t, ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(
                    HartreeEnergies,
                    DipoleMatrix,
                    lasers);

            CrankNicolsonSolver<typename GlobalStateManager::template OperationSpace<1>>
                Solver(
                    Hamiltonian.getMatrix(),
                    TimeMaster::Clock().getTimeStep());

            GlobalStateManager::applyHamiltonian<1>(Solver, m_GlobalStateVector, { { affectedQubit } }, Hamiltonian.getMatrix());
        }
    };

}