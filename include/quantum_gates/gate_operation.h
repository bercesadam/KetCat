#pragma once
#include "gates.h"


namespace KetCat
{
    template<natural_t QubitCount>
    struct GateOperation
    {
        GateType m_type;

        qbit_list_t<QubitCount> m_targets;

        real_t m_theta = 0.0;
    };
}