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

#include "simulation_observer.h"


namespace KetCat
{
    /// @brief Defining these as global constants here, as they work out well and
    /// currently I see no point to expose them ie. in the contructor the the QPU
    /// so it grabs these values directly from here.
    constexpr real_t CrankNicolsonTimeStep = 50; // a.u.
    constexpr natural_t SimuSaveNthFrame = 5E6;

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

		/// @brief Helper for managing the global state vector and mapping between local operations and global state.
		/// As it's a static class, we only need to define the type alias here.
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;

        /// @brief Physical manifold defining the atomic structure and seed states.
        NeutralAtomManifold<Config> m_Manifold;

        /// @brief Full system state vector in the product Hilbert space.
        StateVector<typename GlobalStateManager::FullHilbertSpace> m_GlobalStateVector;

        /// @brief Controller for pulse generation and rotating frame management.
        LaserPulseSequencer<Config, QubitCount> m_laserSequencer;

        /// @brief Helper class which generates the simulation data file output.
        SimulationObserver<QubitCount, Config> m_SimulationObserver;

    public:
        /// @brief Initialize the processor with a specific time-step for TDSE.
        ///
        /// @param simulationOutputFileName Output filename for the simulation data
        constexpr QuantumProcessor(std::string simulationOutputFileName) :
            m_SimulationObserver(m_Manifold, simulationOutputFileName, SimuSaveNthFrame)
        {
            TimeMaster::Clock().init(CrankNicolsonTimeStep);
            StateVector<decltype(m_Manifold)::SingleAtomOperationHilbertSpace>
                OneAtomSeed = m_Manifold.getOperationSeed();
            m_GlobalStateVector = GlobalStateManager::productStateFromSeed(OneAtomSeed);
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
        /// @brief Translate a logical gate into physical hardware instructions.
        ///
        /// @details
        ///    Invokes the GateCompiler to produce a sequence of PhysicalInstructions.
        ///    A single logical gate (e.g., Hadamard) may result in multiple 
        ///    physical pulses or frame rotations.
        template<typename GateOp>
        void executeGate(const GateOp& gate)
        {
			m_SimulationObserver.setSimulationStepName(gateNameToString(gate.m_type) +
                (gate.m_theta > 0.0 ? " (" + std::to_string(gate.m_theta) + ")" : ""));

            std::cout << "Compiling gate: " << gateNameToString(gate.m_type) << ", Theta: " << gate.m_theta << std::endl;

            GateCompiler Compiler;
            auto [PhysicalInstructions, InstructionCount] = Compiler.compile(gate);

            std::cout << "Target qubits: " << gate.m_targets[0] << (InstructionCount > 1 ? ", " + std::to_string(gate.m_targets[1]) : "") << std::endl;
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
        void executeInstruction(const PhysicalInstruction& instruction)
        {
			std::cout << "Executing instruction of type: " << static_cast<int>(instruction.m_type) << std::endl;

            auto PulseEnvelope = m_laserSequencer.calculateLaserEnvelope(instruction);

            // In case of Virtual Z pulses, frame changes are handled internally 
            // by the laser controller and no numerical propagation is required.
            if (!PulseEnvelope)
            {
                return;
            }

            TwoPhotonLaserEnvelope& Envelope = *PulseEnvelope;
            real_t TransitionTimeLimit = Envelope.getTransitionTimeLimit();
            real_t TimeShift = Envelope.getStartTime();

			std::cout << "Starting pulse evolution for instruction. Transition time limit: " << TransitionTimeLimit << " a.u." << std::endl;
			std::cout << "Theta: " << instruction.m_theta << " radians, Phase: " << instruction.m_phase << " radians" << std::endl;

			LaserPulse Pump, Stokes;

            while (TimeMaster::Clock().getCurrentInstructionTime() < TransitionTimeLimit)
            {
                // Obtain current Rabi amplitudes for Pump (Ωp) and Stokes (Ωs)
                std::tie(Pump, Stokes) =
                    Envelope(TimeMaster::Clock().getCurrentInstructionTime() + TimeShift);

                // Map laser fields to the corresponding energy levels in the operation space
                typename MultiRwaRabiHamiltonian<ConfigType::LevelCount>::template laser_array_t Lasers;
                Lasers[ConfigType::Logical0Level] = Pump;
                Lasers[ConfigType::Logical0Level + 1] = Stokes;

                // Propagate the global wavefunction by one time step Δt
                evolveGlobalState(Lasers, instruction.m_targets[0]);

				m_SimulationObserver.exportStep(m_GlobalStateVector, Pump, Stokes);

                TimeMaster::Clock().tick();
            }

			// Ensure that the final state at the end of the pulse is captured
            m_SimulationObserver.exportStep(m_GlobalStateVector, Pump, Stokes, KEYFRAME);

            /// Reset instruction-local timing state for the next pulse
            TimeMaster::Clock().resetCurrentInstructionClock();

			std::cout << "Completed pulse evolution for instruction." << std::endl;
            std::cout << "------------------------------------------" << std::endl;
        }

        /// @brief Perform a single step of the Crank-Nicolson evolution.
        ///
        /// @details
        ///    1. Constructs the local RWA Hamiltonian Ĥ(t).
        ///    2. Builds the unitary propagator U(Δt) via Crank-Nicolson.
        ///    3. Applies U(Δt) to the target qubit in the global state vector.
        void evolveGlobalState(const MultiRwaRabiHamiltonian<ConfigType::LevelCount>::laser_array_t lasers, const natural_t affectedQubit)
        {
            static const std::array<real_t, ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(HartreeEnergies, DipoleMatrix, lasers);

            CrankNicolsonSolver<typename GlobalStateManager::template OperationSpace<1>>
                Solver(Hamiltonian.getMatrix(), TimeMaster::Clock().getTimeStep());

            // Map the local 1-qubit Hamiltonian operation to the global N-qubit state vector
            std::array<natural_t, 1> targets = { affectedQubit };
            GlobalStateManager::applyHamiltonian<1>(Solver, m_GlobalStateVector, targets, Hamiltonian.getMatrix());
        }
    };
}