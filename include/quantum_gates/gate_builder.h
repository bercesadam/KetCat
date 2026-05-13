#pragma once
#include "gate_operation.h"


namespace KetCat
{
    /// @brief Semantic abstraction for quantum gates.
    ///
    /// @details
    ///    Provides a fluent interface for defining gates within a circuit. 
    ///    This class acts as a factory for GateOperation, bridging the gap 
    ///    between abstract gate definitions and their logical application 
    ///    to specific qubit indices (bits).
    ///
    /// @tparam QubitCount The number of qubits the gate acts upon (arity).
    /// @tparam Gate The specific logical GateType (e.g., X, H, RX).
    template<natural_t QubitCount, GateType Gate>
    class QuantumGate
    {
        /// @brief Rotation angle θ in radians.
        /// @details used primarily for parameterized gates like RX(θ), RY(θ), and RZ(θ).
        real_t m_theta = 0.0;

    public:
        constexpr QuantumGate() = default;

        /// @brief Assign a rotation angle to the gate.
        /// 
        /// @param theta Angle in radians.
        /// @return Reference to the gate for method chaining.
        constexpr QuantumGate& withTheta(real_t theta)
        {
            m_theta = theta;
            return *this;
        }

        /// @brief Bind the gate to specific qubit indices.
        ///
        /// @tparam Bits Variadic list of qubit indices.
        /// @param bits The specific indices the gate should target.
        ///
        /// @return A GateOperation POD containing the logic and target indices.
        ///
        /// @details
        ///    Finalizes the logical layer by mapping the abstract gate 
        ///    and its parameters (θ) to physical qubit addresses.
        template<typename... Bits>
        constexpr auto toBits(Bits... bits) const
        {
            static_assert(sizeof...(Bits) == QubitCount,
                "Number of provided bits must match the gate arity.");

            return GateOperation<QubitCount>
            {
                Gate,
                { static_cast<natural_t>(bits)... },
                    m_theta
            };
        }
    };
}