#pragma once
#include "core_types.h"


namespace KetCat
{
    enum class PhysicalInstructionType
    {
        RamanRotation,
        VirtualZ,
        RydbergBlockade
    };

    struct PhysicalInstruction
    {
        PhysicalInstructionType m_type;

        qbit_list_t<2>  m_targets{};

        natural_t m_targetCount = 0;

        real_t m_theta = 0.0;

        real_t m_phase = 0.0;
    };
}