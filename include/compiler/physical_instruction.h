#pragma once
#include <ostream>
#include "core_types.h"


namespace KetCat
{
	/// @brief Enumeration of physical instruction types for neutral atom control.
    enum class PhysicalInstructionType
    {
        // NEUTRAL ATOM INSTRUCTIONS //////////////////////////////////////////////////////////////////////////////////
        
		// @brief Single-qubit Raman rotation between logical states |0⟩ and |1⟩ via an intermediate state.
        RamanRotation,

		// @brief Frame update for a virtual Z rotation, which does not require physical pulses.
        VirtualZ,

		// @brief Two-qubit Rydberg blockade interaction, which applies a conditional phase shift
        // based on the presence of a Rydberg excitation.
        RydbergExcitation,

		// @brief Free evolution under the influence of the atomic Hamiltonian.
        FreeEvolution,

		// TRAPPED ION INSTRUCTIONS ///////////////////////////////////////////////////////////////////////////////////

		// @brief Single-qubit rotation in trapped ion systems
        CarrierRotation,

		// @brief Two-qubit entangling gate in trapped ion systems
        MoelmerSoerensen,
    };


	/// @brief Utility function to convert PhysicalInstructionType enum values to human-readable strings.
    std::string instructionNameToString(PhysicalInstructionType instructionType)
    {
        switch (instructionType)
        {
        case PhysicalInstructionType::RamanRotation: return "RamanRotation";
        case PhysicalInstructionType::VirtualZ: return "VirtualZ";
        case PhysicalInstructionType::RydbergExcitation: return "RydbergExcitation";
        case PhysicalInstructionType::FreeEvolution: return "FreeEvolution";
        default: return "UnknownInstruction";
        }
    }

    /// @brief Placeholder value for the targets array of PhysicalInstruction to be able to determine
	/// the number of targets without the full struct as well, allowing more elegant handling for debugging.
    static constexpr natural_t TARGET_INACTIVE = std::numeric_limits<natural_t>::max();

    struct TimeDependentPhase
    {
        real_t m_time  = -1.0;
        real_t m_phase =  0.0;
    };

    template<natural_t Count>
    using phase_timings_t = std::array< TimeDependentPhase, Count>;

    static constexpr natural_t MAX_PHASE_TIMINGS = 2U;

	/// @brief Concrete instance of a physical control instruction, ready for pulse generation and execution.
    struct PhysicalInstruction
    {
		/// @brief Type of the physical instruction (e.g., RamanRotation, VirtualZ).
        PhysicalInstructionType m_type;

		//// @brief Target qubit indices for the instruction.
        qbit_list_t<2>  m_targets{ TARGET_INACTIVE , TARGET_INACTIVE };

		// @brief Number of valid target indices in m_targets.
        natural_t m_targetCount = 0;

		// @brief Rotation angle θ for Raman rotations or frame updates, in radians.
        real_t m_theta = 0.0;

		// @brief Phase φ for Raman rotations, in radians.
        phase_timings_t<MAX_PHASE_TIMINGS> m_phases;

		friend std::ostream& operator<<(std::ostream& os, const PhysicalInstruction& instruction)
		{
			os << "PhysicalInstruction(Type: " << instructionNameToString(instruction.m_type)
				<< ", Target(s): [";
			for (natural_t i = 0; i < instruction.m_targetCount; ++i)
			{
				os << instruction.m_targets[i];
				if (i < instruction.m_targetCount - 1)
				{
					os << ", ";
				}
			}
			os << "], Theta: " << instruction.m_theta
				<< ", Phase: " << instruction.m_phases[0].m_phase
				<< ")";
			return os;
		}
    };
}