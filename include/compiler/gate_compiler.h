#pragma once
#include <utility>
#include "quantum_gates/gate_operation.h"
#include "physical_instruction.h"


namespace KetCat
{
    /// @brief Maximum number of physical primitives a single logical gate can decompose into.
    constexpr natural_t MaxInstructionsPerGate = 4;

    /// @brief Fixed-size container for compiled physical instructions.
    using physical_instructions_list_t = std::array<PhysicalInstruction, MaxInstructionsPerGate>;

    /// @brief Logical-to-Physical translation layer (Instruction Selection).
    ///
    /// @details
    ///    The GateCompiler maps abstract quantum gates to platform-specific 
    ///    control primitives. It handles decomposition logic such as:
    ///
    ///      • Native Mapping: Direct translation (e.g., RX → RamanRotation).
    ///      • Virtualization: Mapping logical Z-axis rotations to software phase updates.
    ///      • Decomposition: Breaking down composite gates (e.g., H = Rz(π) -> Ry(π/2)).
    class GateCompiler
    {
        natural_t m_UsedInstructionsCount;
        physical_instructions_list_t m_Instructions;

        /// @brief Internal helper to push a new physical primitive onto the local instruction buffer.
        constexpr void append(PhysicalInstructionType type,
            std::initializer_list<natural_t> targets,
            natural_t targetCount,
            real_t theta,
            real_t phase)
        {
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
                m_Instructions[m_UsedInstructionsCount].m_targets[i++] = Target;
            }

            m_UsedInstructionsCount++;
        }

    public:
        /// @brief Compiles a logical gate operation into an ordered sequence of physical instructions.
        ///
        /// @tparam QubitCount Total qubits in the register.
        /// @param op The logical gate operation to be decomposed.
        ///
        /// @return A pair containing the fixed-size instruction array and the actual count of used elements.
        template<natural_t QubitCount>
        constexpr auto compile(const GateOperation<QubitCount>& op)
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
                // Hadamard decomposition: H = Rz(π)Ry(π/2)
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