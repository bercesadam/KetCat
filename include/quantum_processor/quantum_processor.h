#pragma once
#include <array>
#include <variant>
#include <tuple>
#include <utility>
#include <algorithm> 
#include <ranges>
#include <type_traits>

#include "logo.h"

#include "local_space/neutral_atom_manifold.h"
#include "global_space/subspace_operations.h"

#include "quantum_circuit.h"
#include "compiler/gate_compiler.h"

#include "laser/pulse_sequencer.h"

#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "hamiltonian/two_atom_rydberg.h"

#include "solvers/crank_nicolson_solver.h"

#include "simulation_observer.h"


namespace KetCat
{
    /// @brief Defining these as global constants here, as they work out well and
    /// currently I see no point to expose them ie. in the contructor the the QPU
    /// so it grabs these values directly from here.
    constexpr natural_t SimuSaveNthFrame = 1E6;
    constexpr real_t TimeStepsPerInstruction = 1E8;

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
		/// @param initialBitString Bitstring representing the initial state of the qubits
		/// @details This is a private constructor allowing only diagnostic routines to initialize the processor with a specific state.
		/// The public constructor initializes all qubits to |0⟩ by default, corresponding to the ground state of the atoms by default.
        QuantumProcessor(std::string simulationOutputFileName, std::bitset<QubitCount> initialBitString) :
            m_SimulationObserver(m_Manifold, simulationOutputFileName, SimuSaveNthFrame)
        {
            m_GlobalStateVector =
                GlobalStateManager::basisStateFromBitstring(initialBitString,
                    ConfigType::Logical0Level, ConfigType::Logical1Level);

            std::ranges::for_each(std::views::iota(0U, QubitCount) |
                std::views::filter([&](auto i) { return initialBitString.test(i); }),
                [&](auto i) { m_laserSequencer.initializeAtomAsLogical1(i); });
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
			m_SimulationObserver.setSimulationStepName(gateNameToString(gate.m_type) +
                (gate.m_theta > 0.0 ? " (" + std::to_string(gate.m_theta) + ")" : ""));

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
			std::cout << "Executing instruction of type: " << static_cast<int>(instruction.m_type) << std::endl;

            auto PulseEnvelope = m_laserSequencer.calculateLaserEnvelope(instruction);

            // In case of Virtual Z pulses, frame changes are handled internally 
            // by the laser controller and no numerical propagation is required.
            if (!PulseEnvelope)
            {
                return;
            }

            const TwoPhotonLaserEnvelope& Envelope = *PulseEnvelope;
            const real_t TransitionTimeLimit = Envelope.getTransitionTimeLimit();
            const real_t TimeShift = Envelope.getStartTime();
            TimeMaster::Clock().setTimeStep(TransitionTimeLimit / TimeStepsPerInstruction);

			std::cout << "Starting pulse evolution for instruction. Transition time limit: " << TransitionTimeLimit << " a.u. (" <<
                (Units::AtomicTimeToSeconds * TransitionTimeLimit) * 1E9 << " ns" << std::endl;
			std::cout << "Time step: " << TimeMaster::Clock().getTimeStep() << " a.u. (" <<
                (Units::AtomicTimeToSeconds * TimeMaster::Clock().getTimeStep()) * 1E9 << " ns)" << std::endl;
			std::cout << "Theta: " << instruction.m_theta << " radians, Phase: " << instruction.m_phase << " radians" << std::endl;

			LaserPulse Pump, Stokes;

            while (TimeMaster::Clock().getCurrentInstructionTime() < TransitionTimeLimit)
            {
                // Obtain current Rabi amplitudes for Pump (Ωp) and Stokes (Ωs)
                std::tie(Pump, Stokes) =
                    Envelope(TimeMaster::Clock().getCurrentInstructionTime() + TimeShift);

                // Map laser fields to the corresponding energy levels in the operation space
                const natural_t GroundLevelIndex = Envelope.getGroundLevelIndex();
                typename MultiRwaRabiHamiltonian<ConfigType::LevelCount>::template laser_array_t Lasers;
                Lasers[GroundLevelIndex] = Pump;
                Lasers[GroundLevelIndex + 1] = Stokes;

                // Propagate the global wavefunction by one time step Δt
				// Depending on the instruction type, we evolve either a single qubit (Raman rotation) or two qubits (Rydberg blockade).
                if (instruction.m_type == PhysicalInstructionType::RamanRotation)
                {
                    evolveOneQubitGlobalState(Lasers, instruction.m_targets[0]);
                }
                else if (instruction.m_type == PhysicalInstructionType::RydbergExcitation)
                {
                    if constexpr (QubitCount >= 2)
                    {
                        //evolveOneQubitGlobalState(Lasers, instruction.m_targets[0]);
                        evolveTwoQubitGlobalState(Lasers, instruction.m_targets[0], instruction.m_targets[1]);
                    }
                }
                else { }
               
				/// Capture the current state and laser configuration for visualization/export.
				m_SimulationObserver.exportStep(m_GlobalStateVector,
                    instruction.m_targets, Pump, Stokes);

                TimeMaster::Clock().tick();
            }

			// Ensure that the final state at the end of the pulse is captured
            m_SimulationObserver.exportStep(m_GlobalStateVector,
                instruction.m_targets, Pump, Stokes, KEYFRAME);

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
        void evolveOneQubitGlobalState(
            const MultiRwaRabiHamiltonian<ConfigType::LevelCount>::laser_array_t lasers,
            const natural_t affectedQubit)
        {
            static const std::array<real_t, ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(HartreeEnergies, DipoleMatrix, lasers);

            static CrankNicolsonSolver<ConfigType::LevelCount, LinearSolverBackend::ThomasTridiagonal> Solver;

			Hamiltonian.updateMainDiagonal(lasers);
			Hamiltonian.updateOffDiagonal(lasers);
			Solver.updateMatrices(Hamiltonian.getMatrix(), TimeMaster::Clock().getTimeStep());

            // Map the local 1-qubit Hamiltonian operation to the global N-qubit state vector
            std::array<natural_t, 1> targets = { affectedQubit };
             GlobalStateManager::template performTimeEvolution<1>(Solver, m_GlobalStateVector, targets);
        }

        void evolveTwoQubitGlobalState(
            const MultiRwaRabiHamiltonian<ConfigType::LevelCount>::laser_array_t lasers,
            const natural_t controlAtom, const natural_t targetAtom)
            requires (QubitCount >= 2)
        {
            static const std::array<real_t, ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const square_matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            static MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                SingleAtomExcitation(HartreeEnergies, DipoleMatrix, lasers);

            static TwoAtomRydbergBlockade<ConfigType::LevelCount>
                RydbergBlockade(Units::MeterToAtomicLength * 1E-9,
                    ConfigType::RydbergLevel,
                    HartreeEnergies,
				DipoleMatrix);

            // Use FiveBandGaussianElimination backend for two-qubit dense Hamiltonians
            static CrankNicolsonSolver<ConfigType::LevelCount, LinearSolverBackend::FiveBandGaussianElimination> Solver;

            SingleAtomExcitation.updateMainDiagonal(lasers);
            SingleAtomExcitation.updateOffDiagonal(lasers);
            auto SingleAtomHamiltonian = SingleAtomExcitation.getMatrix();

            RydbergBlockade.updateMatrix(SingleAtomHamiltonian);

            Solver.updateMatrices(RydbergBlockade.getMatrix(), TimeMaster::Clock().getTimeStep());

            // Map the local 1-qubit Hamiltonian operation to the global N-qubit state vector
            std::array<natural_t, 2> targets = { controlAtom, targetAtom };
            GlobalStateManager::template performTimeEvolution<2>(Solver, m_GlobalStateVector, targets);
        }

    private:
       // Grant access to the internal state and methods of QuantumProcessor for diagnostic purposes.
        template<natural_t, NeutralAtomTypeConfig>
        friend class QPUDiagnostics;
    };
}