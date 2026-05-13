#pragma once
#include "gates.h"


namespace KetCat
{
	/// @brief Concrete instance of a quantum gate operation, ready for compilation and execution.
    template<natural_t QubitCount>
    struct GateOperation
    {
        GateType m_type;

        qbit_list_t<QubitCount> m_targets;

        real_t m_theta = 0.0;
    };
}