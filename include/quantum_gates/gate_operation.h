#pragma once
#include "gates.h"


namespace KetCat
{
	/// @brief Concrete instance of a quantum gate operation, ready for compilation and execution.
    template<natural_t QubitCount>
    struct GateOperation
    {
		// @brief Type of the logical gate (e.g., X, H, RX).
        GateType m_type;

		/// @brief List of target qubit indices for the gate operation.
        qbit_list_t<QubitCount> m_targets;

		/// @brief Rotation angle θ for parameterized gates, in radians.
        real_t m_theta = 0.0;
    };
}