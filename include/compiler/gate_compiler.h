#pragma once
#include <array>
#include <utility>
#include "quantum_gates/gate_operation.h"
#include "physical_instruction.h"

namespace KetCat
{
    /// @brief Maximum number of physical primitives a single logical gate can decompose into.
    constexpr natural_t MaxInstructionsPerGate = 10;

    /// @brief Fixed-size container for compiled physical instructions.
    using physical_instructions_list_t =
        std::array<PhysicalInstruction, MaxInstructionsPerGate>;

    /// @brief Logical-to-Physical translation layer (Instruction Selection).
    ///
    /// @details
    ///     Maps logical quantum gates to native hardware primitives.
    ///
    ///     Supports:
    ///         • Native gate mappings
    ///         • Virtual Z rotations
    ///         • Gate decompositions
    ///         • Compile-time lowering paths
    class GateCompiler
    {
    public:
        using GateCompilerResult =
            std::pair<physical_instructions_list_t, natural_t>;

    private:
        natural_t m_UsedInstructionsCount = 0;
        physical_instructions_list_t m_Instructions = {};

        /// @brief Append a native physical instruction.
        constexpr void append(
            PhysicalInstructionType type,
            std::initializer_list<natural_t> targets,
            natural_t targetCount,
            real_t theta,
            real_t phase)
        {
            PhysicalInstruction& Instruction =
                m_Instructions[m_UsedInstructionsCount];

            Instruction =
            {
                type,
                { TARGET_INACTIVE, TARGET_INACTIVE },
                targetCount,
                theta,
                phase
            };

            natural_t i = 0;

            for (natural_t Target : targets)
            {
                Instruction.m_targets[i++] = Target;
            }

            ++m_UsedInstructionsCount;
        }

        /// @brief Append all instructions from another compiler result.
        constexpr void append(const GateCompilerResult& result)
        {
            for (natural_t i = 0; i < result.second; ++i)
            {
                m_Instructions[m_UsedInstructionsCount++] =
                    result.first[i];
            }
        }

    private:

        /// @brief Compile-time gate lowering path.
        template<GateType Type, natural_t QubitCount>
        constexpr void compileGate(const GateOperation<QubitCount>& op)
        {
            ///////////////////////////////////////////
            // -- BLOCH ROTATIONS --
            //
            // RX
            //
            if constexpr (Type == GateType::RX)
            {
                append(
                    PhysicalInstructionType::RamanRotation,
                    { op.m_targets[0] },
                    1,
                    op.m_theta,
                    0.0);
            }

            //
            // RY
            //
            else if constexpr (Type == GateType::RY)
            {
                append(
                    PhysicalInstructionType::RamanRotation,
                    { op.m_targets[0] },
                    1,
                    op.m_theta,
                    -ConstexprMath::Pi / 2.0);
                    }

            //
            // RZ
            //
            else if constexpr (Type == GateType::RZ)
            {
                append(
                    PhysicalInstructionType::VirtualZ,
                    { op.m_targets[0] },
                    1,
                    op.m_theta,
                    0.0);
            }
            

            ///////////////////////////////////////////
            // -- PAULI GATES (“π” rotations) --

            // Pauli X
            else if constexpr (Type == GateType::X)
            {
                const natural_t target = op.m_targets[0];
                compileGate<GateType::RX>(GateOperation<1>{
                    GateType::RX,
                    { target },
                        ConstexprMath::Pi
                });
            }

            // Pauli Y
            else if constexpr (Type == GateType::Y)
            {
                const natural_t target = op.m_targets[0];
                compileGate<GateType::RY>(GateOperation<1>{
                    GateType::RY,
                    { target },
                        ConstexprMath::Pi
                });
            }

            // Pauli Z
            else if constexpr (Type == GateType::Z)
            {
                const natural_t target = op.m_targets[0];
                compileGate<GateType::RZ>(GateOperation<1>{
                    GateType::RZ,
                    { target },
                        ConstexprMath::Pi
                });
            }

            ///////////////////////////////////////////
            //
            // Hadamard
            //
            else if constexpr (Type == GateType::H)
            {
                // H = Z * RY(pi/2)

                const natural_t target = op.m_targets[0];

                compileGate<GateType::Z>(GateOperation<1>{
                    GateType::Z,
                    { target }
                });

                compileGate<GateType::RY>(GateOperation<1>{
                    GateType::RY,
                    { target },
                        ConstexprMath::Pi / 2.0
                });
            }

            ///////////////////////////////////////////
            //
            // Controlled-Z (CPHASE / CZ)
            //
            else if constexpr (Type == GateType::CZ)
            {
                //
                // Native Rydberg blockade controlled phase.
                //
                // Applies:
                //
                //      |11> -> -|11>
                //
                // while leaving all other computational basis
                // states unchanged.
                //
                const natural_t ctrl = op.m_targets[0];
                const natural_t trgt = op.m_targets[1];

                constexpr real_t phaseShift = ConstexprMath::Pi + (ConstexprMath::Pi / 2.23606797749979); // 2.236067... = sqrt(5)

                append(
                    PhysicalInstructionType::RydbergExcitation,
                    { ctrl, trgt },
                    2,
                    ConstexprMath::Pi,
                    0.0
                );

                append(
                    PhysicalInstructionType::RydbergExcitation,
                    { ctrl, trgt },
                    2,
                    ConstexprMath::Pi,
                    phaseShift
                );
            }

            ///////////////////////////////////////////
            //
            // Controlled-X (CNOT)
            //
            else if constexpr (Type == GateType::CX)
            {
                //
                // CNOT decomposition:
                //
                //      CX = H(target) · CZ · H(target)
                //
                const natural_t ctrl = op.m_targets[0];
                const natural_t trgt = op.m_targets[1];

                // H(target)
                compileGate<GateType::H>(GateOperation<1>{
                    GateType::H,
                    { trgt }
                });

                // CZ(control, target)
                compileGate<GateType::CZ>(GateOperation<2>{
                    GateType::CZ,
                    { ctrl, trgt }
                });

                // H(target)
                compileGate<GateType::H>(GateOperation<1>{
                    GateType::H,
                    { trgt }
                });
            }

            //
            // Controlled-Controlled-X (CCX / Toffoli)
            //
            else if constexpr (Type == GateType::CCX)
            {
                const natural_t control1 = op.m_targets[0];
                const natural_t control2 = op.m_targets[1];
                const natural_t target = op.m_targets[2];

                //
                // CCX decomposition:
                //
                // CCX = H(t) · CCZ · H(t)
                //

                //
                // H(target)
                //
                compileGate<GateType::H>(
                    GateOperation<QubitCount>
                {
                    GateType::H,
                    { target }
                });

                //
                // CCZ decomposition:
                // CCZ = CZ(c1,t) · T(t) · CZ(c2,t) · T†(t) · CZ(c1,t)
                //
                // (standard Clifford+T decomposition without ancilla)
                //

                compileGate<GateType::CZ>(
                    GateOperation<2>
                {
                    GateType::CZ,
                    { control1, target }
                });

                compileGate<GateType::RZ>(
                    GateOperation<QubitCount>
                {
                    GateType::RZ,
                    { target },
					ConstexprMath::Pi / 4.0
                });

                compileGate<GateType::CZ>(
                    GateOperation<2>
                {
                    GateType::CZ,
                    { control2, target }
                });

                compileGate<GateType::RZ>(
                    GateOperation<QubitCount>
                {
                    GateType::RZ,
                    { target },
                    -ConstexprMath::Pi / 4.0
                });

                compileGate<GateType::CZ>(
                    GateOperation<2>
                {
                    GateType::CZ,
                    { control1, target }
                });

                //
                // H(target)
                //
                compileGate<GateType::H>(
                    GateOperation<QubitCount>
                {
                    GateType::H,
                    { target }
                });
            }
        }

    public:

        /// @brief Compile a logical gate into physical primitives.
        template<natural_t QubitCount>
        constexpr GateCompilerResult compile(
            const GateOperation<QubitCount>& op)
        {
            m_UsedInstructionsCount = 0;
            m_Instructions = {};

            switch (op.m_type)
            {
            case GateType::X:
                compileGate<GateType::X>(op);
                break;

            case GateType::Y:
                compileGate<GateType::Y>(op);
                break;

            case GateType::Z:
                compileGate<GateType::Z>(op);
                break;

            case GateType::RX:
                compileGate<GateType::RX>(op);
                break;

            case GateType::RY:
                compileGate<GateType::RY>(op);
                break;

            case GateType::RZ:
                compileGate<GateType::RZ>(op);
                break;

            case GateType::H:
                compileGate<GateType::H>(op);
                break;

            case GateType::CZ:
                compileGate<GateType::CZ>(op);
                break;

            case GateType::CX:
                compileGate<GateType::CX>(op);
                break;

            case GateType::CCX:
                compileGate<GateType::CCX>(op);
				break;

            default:
                break;
            }

            return
            {
                m_Instructions,
                m_UsedInstructionsCount
            };
        }
    };
}