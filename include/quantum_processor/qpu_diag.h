#pragma once
#include "quantum_processor.h"


namespace KetCat
{
	/// @brief Diagnostic wrapper for auditing and driving physical simulation states of a neutral-atom QPU.
	/// @tparam QubitCount The total number of physical atoms layouted in the architecture.
	/// @tparam Config Hardware topology and specific atomic species configuration traits.
	template <natural_t QubitCount, NeutralAtomTypeConfig Config>
	class QPUDiagnostics
	{
		using QPUType = QuantumProcessor<QubitCount, Config>;

		/// @brief Underlyling active hardware simulator instance.
		QPUType m_processor;

		/// @brief Construct a diagnostic session tracking a specific QPU register setup.
		/// @param simulationOutputFileName Path to serialize detailed execution and execution waveforms.
		/// @param initialBitString Binary layout mapping the computational basis state |b₁b₂...bₙ⟩.
		QPUDiagnostics(std::string simulationOutputFileName, std::bitset<QubitCount> initialBitString)
			: m_processor(simulationOutputFileName, initialBitString)
		{
			std::cout << "KETCAT AB INITIO NEUTRAL ATOM QUANTUM COMPUTER SIMULATOR v3.0 - Diagnostic Session initiated." << std::endl;
		}

	public:
		/// @brief Factory method creating a QPU Diagnostic instance pre-initialized to a computational state.
		/// @param simulationOutputFileName Path to serialize detailed execution and execution waveforms.
		/// @param initialBitString Binary layout mapping the computational basis state |b₁b₂...bₙ⟩.
		/// @return Fully populated diagnostics manager instance.
		static QPUDiagnostics createQPUWithInitialState(std::string simulationOutputFileName, std::bitset<QubitCount> initialBitString)
		{
			return QPUDiagnostics(simulationOutputFileName, initialBitString);
		}

		/// @brief Exposes a mutable reference to the active internal quantum processor simulator.
		/// @return Reference to the QuantumProcessor core engine.
		QPUType& QPU()
		{
			return m_processor;
		}

		/// @brief Returns the full uncompressed state vector representing the multi-atom Hilbert space.
		const auto& getGlobalStateVector()
		{
			return m_processor.m_GlobalStateManager.getStateVector();
		}

		/// @brief Directly dispatches a native assembly-level physical instruction onto the neutral-atom array.
		void executePhysicalInstruction(const PhysicalInstruction& instruction)
		{
			m_processor.executeInstruction(instruction);
		}
	};
}