#pragma once
#include "core_types.h"


namespace KetCat
{
	/// @brief Enumeration of physical instruction types for neutral atom control.
    enum class PhysicalInstructionType
    {
        RamanRotation,
        VirtualZ,
        RydbergBlockade
    };

	/// @brief Concrete instance of a physical control instruction, ready for pulse generation and execution.
    struct PhysicalInstruction
    {
		/// @brief Type of the physical instruction (e.g., RamanRotation, VirtualZ).
        PhysicalInstructionType m_type;

		//// @brief Target qubit indices for the instruction.
        qbit_list_t<2>  m_targets{};

		// @brief Number of valid target indices in m_targets.
        natural_t m_targetCount = 0;

		// @brief Rotation angle θ for Raman rotations or frame updates, in radians.
        real_t m_theta = 0.0;

		// @brief Phase φ for Raman rotations, in radians.
        real_t m_phase = 0.0;
    };
}