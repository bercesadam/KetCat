#pragma once
#include "quantum_processor.h"


namespace KetCat
{
	template <natural_t QubitCount, NeutralAtomTypeConfig Config>
	class QPUDiagnostics
	{
		using QPUType = QuantumProcessor<QubitCount, Config>;

		QPUType m_processor;

		QPUDiagnostics(std::string simulationOutputFileName, std::bitset<QubitCount> initialBitString)
			: m_processor(simulationOutputFileName, initialBitString)
		{
			std::cout << "KETCAT AB INITIO NEUTRAL ATOM QUANTUM COMPUTER SIMULATOR v3.0 - Diagnostic Session initiated." << std::endl;
		}

	public:
		static QPUDiagnostics createQPUWithInitialState(std::string simulationOutputFileName, std::bitset<QubitCount> initialBitString)
		{
			return QPUDiagnostics(simulationOutputFileName, initialBitString);
		}

		QPUType& QPU()
		{
			return m_processor;
		}

		const auto& getGlobalStateVector()
		{
			return m_processor.m_GlobalStateManager.getStateVector();
		}

		void executePhysicalInstruction(const PhysicalInstruction& instruction)
		{
			m_processor.executeInstruction(instruction);
		}
	};
}