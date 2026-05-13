#pragma once
#include "core_types.h"
#include "quantum_gates/gate_builder.h"

namespace KetCat
{
    using OperationVariant = std::variant<GateOperation<1>, GateOperation<2>>;

    template<natural_t QubitCount, natural_t OpCount = 0>
    class QuantumCircuit
    {
        template <natural_t Count>
        using operation_array_t = std::array<OperationVariant, Count>;

        operation_array_t<OpCount> m_operations;

    public:
        constexpr QuantumCircuit() = default;

        constexpr QuantumCircuit(const std::array<OperationVariant, OpCount>& ops)
            : m_operations(ops)
        {}

        constexpr const auto& operations() const noexcept
        {
            return m_operations;
        }

        template<typename GateOp>
        constexpr auto withGate(const GateOp& op) const
        {
            std::array<OperationVariant, OpCount + 1> UpdatedArray;
            std::copy_n(m_operations.begin(), OpCount, UpdatedArray.begin());

            UpdatedArray[OpCount] = op;

            return QuantumCircuit<QubitCount, OpCount + 1>(UpdatedArray);
        }

        template<typename... Ops>
        constexpr auto withGates(Ops... ops) const
        {
            auto result = *this;

            ((result = result.withGate(ops)), ...);

            return result;
        }

        template<typename GateTag, natural_t N>
        constexpr auto applySingleQubitGate(const std::array<natural_t, N>& targets) const
        {
            auto result = *this;

            for (natural_t i = 0; i < N; ++i)
            {
                result = result.withGate(
                    QuantumGate<1, GateTag>()
                    .toBits(targets[i]));
            }

            return result;
        }
    };
}