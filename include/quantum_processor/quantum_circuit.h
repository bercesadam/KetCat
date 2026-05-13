#pragma once
#include <variant>
#include "core_types.h"
#include "quantum_gates/gate_builder.h"

namespace KetCat
{
    /// @brief Type-safe variant representing the supported arity of gate operations.
    using OperationVariant = std::variant<GateOperation<1>, GateOperation<2>>;

    /// @brief Compile-time representation of a quantum circuit.
    ///
    /// @details
    ///    The QuantumCircuit class uses a functional, immutable design to build a 
    ///    sequence of quantum operations at compile-time. Each gate addition 
    ///    results in a new type instance with an incremented GateCount.
    ///
    ///    Design Architecture:
    ///      - Logical Representation: Stores a sequence of GateOperations.
    ///      - Fluent DSL Interface: Supports chaining via withGates and recursion-based 
    ///        template expansion.
    ///
    /// @tparam QubitCount The total width of the quantum register (wires).
    /// @tparam GateCount The current number of logical operations in the circuit depth.
    template<natural_t QubitCount, natural_t GateCount = 0>
    class QuantumCircuit
    {
        template <natural_t Count>
        using operation_array_t = std::array<OperationVariant, Count>;

        /// @brief Sequential storage of operations.
        operation_array_t<GateCount> m_Gates;

    public:
        constexpr QuantumCircuit() = default;

        /// @brief Internal constructor for creating extended circuits from existing gate arrays.
        constexpr QuantumCircuit(const std::array<OperationVariant, GateCount>& gates)
            : m_Gates(gates)
        {}

        /// @brief Access the ordered list of gate operations.
        constexpr const auto& operations() const noexcept
        {
            return m_Gates;
        }

        /// @brief Add multiple gates to the circuit.
        ///
        /// @tparam Gates Variadic list of GateOperation types.
        /// @param gates Instances of gates to append.
        /// @return A new QuantumCircuit instance containing the additional gates.
        template<typename... Gates>
        constexpr auto withGates(Gates... gates) const
        {
            return applyImpl(*this, gates...);
        }

        /// @brief Utility to apply a single-qubit gate across multiple target indices.
        ///
        /// @tparam GateTag The logical GateType to apply.
        /// @tparam N The number of targets provided.
        /// @param targets An array of qubit indices.
        /// @return A new QuantumCircuit instance with N new gates.
        template<typename GateTag, natural_t N>
        constexpr auto applySingleQubitGate(const std::array<natural_t, N>& targets) const
        {
            if constexpr (N == 0)
            {
                return *this;
            }

            return applyTargetsImpl<GateTag>(*this, targets, std::make_index_sequence<N>{});
        }

    private:
        /// @brief Implementation helper to expand target indices into individual gate operations.
        template<typename GateTag, typename Circuit, std::size_t... I, natural_t N>
        static constexpr auto applyTargetsImpl(Circuit c, const std::array<natural_t, N>& targets, std::index_sequence<I...>)
        {
            return applyImpl(c, QuantumGate<1, GateTag>().toBits(targets[I])...);
        }

        /// @brief Base case for recursive gate application.
        template<typename Circuit>
        static constexpr auto applyImpl(Circuit c)
        {
            return c;
        }

        /// @brief Recursive implementation for gate chaining.
        template<typename Circuit, typename Gate, typename... Rest>
        static constexpr auto applyImpl(Circuit c, const Gate& gate, const Rest&... rest)
        {
            return applyImpl(c.addGate(gate), rest...);
        }

        /// @brief Creates a new circuit instance by appending a single gate.
        ///
        /// @details 
        ///    Performs a compile-time copy of the existing gate array and appends 
        ///     the new operation at index [GateCount].
        template<typename Gate>
        constexpr auto addGate(const Gate& gate) const
        {
            std::array<OperationVariant, GateCount + 1> UpdatedArray;
            std::copy_n(m_Gates.begin(), GateCount, UpdatedArray.begin());

            UpdatedArray[GateCount] = gate;

            return QuantumCircuit<QubitCount, GateCount + 1>(UpdatedArray);
        }

		// Granting other QuantumCircuit template instances access to private members for execution.
        template<natural_t, natural_t>
        friend class QuantumCircuit;
    };
}