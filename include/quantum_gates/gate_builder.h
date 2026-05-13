#pragma once
#include "gate_operation.h"


namespace KetCat
{
    template<natural_t QubitCount, GateType Gate>
    class QuantumGate
    {
        /// @brief Parameter for rotation gates (RX, RY, RZ). Ignored for non-parameterized gates.
        real_t m_theta = 0.0;

    public:
        constexpr QuantumGate() = default;

        constexpr QuantumGate& withTheta(real_t theta)
        {
            m_theta = theta;
            return *this;
        }

        template<typename... Bits>
        constexpr auto toBits(Bits... bits) const
        {
            static_assert(sizeof...(Bits) == QubitCount);

            return GateOperation<QubitCount>
            {
                Gate,
                { static_cast<natural_t>(bits)...},
                m_theta
            };
        }
    };
}