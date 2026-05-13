#pragma once
#include <utility>
#include "quantum_gates/gate_operation.h"
#include "physical_instruction.h"


namespace KetCat
{
    constexpr natural_t MaxInstructionsPerGate = 4;
    using physical_instructions_list_t = std::array<PhysicalInstruction, MaxInstructionsPerGate>;

    class GateCompiler
    {
        natural_t m_UsedInstructionsCount;
        physical_instructions_list_t m_Instructions;

        void append(PhysicalInstructionType type,
            std::initializer_list<natural_t> targets,
            natural_t targetCount,
            real_t theta,
            real_t phase)
        {
            m_UsedInstructionsCount++;
            m_Instructions[m_UsedInstructionsCount] =
            {
                type,
                {},
                targetCount,
                theta,
                phase
            };

            natural_t i = 0;
            for (auto Target : targets)
            {
                m_Instructions[m_UsedInstructionsCount - 1].m_targets[i++] = Target;
            }
        }

    public:
        template<natural_t QubitCount>
        static constexpr auto compile(const GateOperation<QubitCount>& op)
        {
            m_UsedInstructionsCount = 0;
            m_Instructions = {};

            switch (op.m_type)
            {
                case GateType::X:
                {
                    append(
                        PhysicalInstructionType::RamanRotation,
                        { op.m_targets[0] },
                        1,
                        ConstexprMath::Pi,
                        0.0);

                    break;
                }

                case GateType::Y:
                {
                    append(
                        PhysicalInstructionType::RamanRotation,
                        { op.m_targets[0] },
                        1,
                        ConstexprMath::Pi,
                        ConstexprMath::Pi / 2.0);

                    break;
                }

                case GateType::RX:
                {
                    append(
                        PhysicalInstructionType::RamanRotation,
                        { op.m_targets[0] },
                        1,
                        op.m_theta,
                        0.0);

                    break;
                }

                case GateType::RY:
                {
                    append(
                        PhysicalInstructionType::RamanRotation,
                        { op.m_targets[0] },
                        1,
                        op.m_theta,
                        ConstexprMath::Pi / 2.0);

                    break;
                }

                case GateType::RZ:
                {
                    append(
                        PhysicalInstructionType::VirtualZ,
                        { op.m_targets[0] },
                        1,
                        op.m_theta,
                        0.0);

                    break;
                }

                case GateType::H:
                {
                    append(
                        PhysicalInstructionType::VirtualZ,
                        { op.m_targets[0] },
                        1,
                        ConstexprMath::Pi,
                        0.0);

                    append(
                        PhysicalInstructionType::RamanRotation,
                        { op.m_targets[0] },
                        1,
                        ConstexprMath::Pi / 2.0,
                        ConstexprMath::Pi / 2.0);

                    break;
                }

                case GateType::CX:
                {
                    append(
                        PhysicalInstructionType::RydbergBlockade,
                        {
                            op.m_targets[0],
                            op.m_targets[1]
                        },
                        2,
                        0.0,
                        0.0);

                    break;
                }

                default:
                    break;
            }

            return std::pair
            {
                m_Instructions,
                m_UsedInstructionsCount
            };
        }
    };
}
