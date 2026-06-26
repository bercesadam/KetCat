#pragma once
#include <variant>
#include <type_traits>

#include "logo.h"
#include "simulation_observer.h"

#include "quantum_circuit.h"
#include "compiler/gate_compiler.h"
#include "laser/pulse_sequencer.h"
#include "global_space/global_state_manager.h"


namespace KetCat
{
    /// @brief Defining these as global constants here, as they work out well and
    /// currently I see no point to expose them ie. in the contructor the the QPU
    /// so it grabs these values directly from here.
    constexpr natural_t SimuSaveNthFrame = 1E6;
    constexpr real_t TimeStepsPerInstruction = 5E7;

	/// @brief Forward declare Diagnostic class for friend declaration.
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class QPUDiagnostics;

    /// @brief Main control logic/orchestraion of the complete neutral atom quantum computer simulation stack.
    ///
    /// @details
	///     The architecture consists of three main, well separated abstraction layers:
    ///      1. Logical Layer:   Abstract gates are received via QuantumCircuit.
    ///      2. Physical Layer:  GateCompiler decomposes gates into PhysicalInstructions.
    ///      3. Execution Layer: LaserPulseSequencer generates envelopes from PhysicalInstructions for TDSE propagation.
    ///
    /// @tparam QubitCount Number of qubits in the simulated register.
    /// @tparam Config Physical atom configuration (levels, dipoles, and logical mapping).
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class QuantumProcessor
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;

        /// @brief Full system state vector in the product Hilbert space.
        GlobalStateManager<QubitCount, Config> m_GlobalStateManager;

        /// @brief Controller for pulse generation and rotating frame management.
        LaserPulseSequencer<Config, QubitCount> m_laserSequencer;

        /// @brief Helper class which generates the simulation data file output.
        SimulationObserver<QubitCount, Config> m_SimulationObserver;

    public:
        /// @brief Initialize the processor with a specific time-step for TDSE.
        ///
        /// @param simulationOutputFileName Output filename for the simulation data
        QuantumProcessor(std::string simulationOutputFileName)
            : QuantumProcessor(simulationOutputFileName, std::bitset<QubitCount>{})
        {  
			printLogo();
        }

        /// @brief Execute a logical quantum circuit.
        ///
        /// @details 
        ///    Iterates through the circuit's variant-based operations, 
        ///    dispatching each to the compiler for physical translation.
        template<natural_t OpCount>
        void execute(const QuantumCircuit<QubitCount, OpCount>& circuit)
        {
            for (const auto& op : circuit.operations())
            {
                std::visit([&](const auto& gate) { executeGate(gate); }, op);
            }
        }

    private:
		/// @brief Internal constructor for initializing the processor with a specific initial state.
		/// @param simulationOutputFileName Output filename for the simulation data
		/// @param initialState Bitstring representing the initial state of the qubits
		/// @details This is a private constructor allowing only diagnostic routines to initialize the processor with a specific state.
		/// The public constructor initializes all qubits to |0⟩ by default, corresponding to the ground state of the atoms by default.
        QuantumProcessor(std::string simulationOutputFileName, std::bitset<QubitCount> initialState) :
            m_GlobalStateManager(initialState),
            m_SimulationObserver(m_GlobalStateManager.getManifold(), simulationOutputFileName, SimuSaveNthFrame)
        {
        }

        /// @brief Translate a logical gate into physical hardware instructions.
        ///
        /// @details
        ///    Invokes the GateCompiler to produce a sequence of PhysicalInstructions.
        ///    A single logical gate (e.g., Hadamard) may result in multiple 
        ///    physical pulses or frame rotations.
        template<typename GateOp>
        void executeGate(const GateOp& gate)
        {
			m_SimulationObserver.appendSimulationStepName(gateNameToString(gate.m_type) +
				(gate.m_theta > 0.0 ? " (" + std::to_string(gate.m_theta) + "); " : "; "));

            std::cout << "Starting gate compilation: " << gateNameToString(gate.m_type) << std::endl;

            GateCompiler Compiler;
            auto [PhysicalInstructions, InstructionCount] = Compiler.compile(gate);

            std::cout << "Target qubits: " << gate.m_targets[0] <<
                (gate.m_targets.size() > 1 ? ", " + std::to_string(gate.m_targets[1]) : "") << std::endl;
			std::cout << "Generated " << InstructionCount << " physical instructions." << std::endl;

            for (natural_t i = 0; i < InstructionCount; ++i)
            {
                executeInstruction(PhysicalInstructions[i]);
            }
        }

        /// @brief Drive the system state based on a physical instruction.
        ///
        /// @details
        ///    High-level logic for handling instruction types:
        ///      • Virtual Instructions: Updates internal sequencer state; returns std::nullopt.
        ///      • Physical Pulses: Evaluates the LaserEnvelope and evolves the state via TDSE.
        ///
        ///    The evolution loop propagates Ψ(t) using the time-dependent Hamiltonian:
        ///      i ∂/∂t |Ψ⟩ = Ĥ(t)|Ψ⟩
        void executeInstruction(PhysicalInstruction& instruction)
        {
            std::cout << "Executing instruction: " << instruction << std::endl;
            m_SimulationObserver.appendSimulationStepName(instructionNameToString(instruction.m_type) + ", θ=" +
                std::to_string(instruction.m_theta) + ", φ=" + std::to_string(instruction.m_phases[0].m_phase) + "; ");

            auto PulseEnvelope = m_laserSequencer.calculateLaserEnvelope(instruction);

            // In case of Virtual Z pulses, frame changes are handled internally 
            // by the laser controller and no numerical propagation is required.
            if (!PulseEnvelope)
            {
                return;
            }

            const TwoPhotonLaserEnvelope& Envelope = *PulseEnvelope;
            const real_t TransitionTimeLimit = Envelope.getTransitionTimeLimit();
            TimeMaster::Clock().setTimeStep(TransitionTimeLimit / TimeStepsPerInstruction);

			std::cout << "Starting pulse evolution for instruction. Transition time limit: " << TransitionTimeLimit << " a.u. (" <<
                (Units::AtomicTimeToSeconds * TransitionTimeLimit) * 1E9 << " ns" << std::endl;
			std::cout << "Time step: " << TimeMaster::Clock().getTimeStep() << " a.u. (" <<
                (Units::AtomicTimeToSeconds * TimeMaster::Clock().getTimeStep()) * 1E9 << " ns)" << std::endl;
			std::cout << "Theta: " << instruction.m_theta << " radians, Phase: " << instruction.m_phases[0].m_phase << " radians" << std::endl;

			TwoPhotonDrive Lasers;

            while (TimeMaster::Clock().getCurrentInstructionTime() < TransitionTimeLimit)
            {
                // Obtain current laser drive parameters
                Lasers = Envelope(TimeMaster::Clock().getCurrentInstructionTime());

                // Propagate the global wavefunction by one time step Δt
				// Depending on the instruction type, we evolve either a single qubit (Raman rotation) or two qubits (Rydberg blockade).
                if (instruction.m_type == PhysicalInstructionType::RamanRotation)
                {
                    m_GlobalStateManager.evolveOneQubitGlobalState(Lasers, instruction.m_targets[0]);
                }
                else if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
                {
                    if constexpr (QubitCount >= 2)
                    {
                        m_GlobalStateManager.evolveTwoQubitGlobalState(Lasers, instruction.m_targets[0], instruction.m_targets[1]);
                    }
                }
                else { }
               
				/// Capture the current state and laser configuration for visualization/export.
				m_SimulationObserver.exportStep(m_GlobalStateManager.getStateVector(),
                    instruction.m_targets, Lasers.m_pump, Lasers.m_stokes);

                TimeMaster::Clock().tick();
            }

			// Ensure that the final state at the end of the pulse is captured
            m_SimulationObserver.exportStep(m_GlobalStateManager.getStateVector(),
                instruction.m_targets, Lasers.m_pump, Lasers.m_stokes, KEYFRAME);

            /// Reset instruction-local timing state for the next pulse
            TimeMaster::Clock().resetCurrentInstructionClock();

            // Reset simulation step title
            m_SimulationObserver.resetSimulationStepName();

			std::cout << "Completed pulse evolution for instruction." << std::endl;
            std::cout << "------------------------------------------" << std::endl;
        }

    private:
       // Grant access to the internal state and methods of QuantumProcessor for diagnostic purposes.
        template<natural_t, NeutralAtomTypeConfig>
        friend class QPUDiagnostics;
    };
}