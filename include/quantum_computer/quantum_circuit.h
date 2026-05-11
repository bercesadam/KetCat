#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @file
    /// @brief Small executor and convenience API for composing and running quantum gates.

	/// @brief Forward declaration of QuantumCircuit for friend declaration.
    template<natural_t QBitCount>
    class QuantumCircuit;

    /// @brief Concept for types that behave like a gate: they are callable with a state vector.
    /// @tparam GateType  Type to test.
    ///
    /// The requirement checks one instantiation (with StateCount == 1) and uses that as a
    /// proxy for gate-like behaviour. Implementers should ensure their operator() is
    /// templated or overloaded appropriately for various `StateCount` values.
    template<typename GateType>
    concept QuantumGateLike =
        requires(GateType g)
    {
        // Try invoking the gate with a sample 1-qubit state vector; the returned type
        // must be a state_vector_t<1>. This models the "callable that maps state->state".
        { g(StateVector<FiniteHilbertSpace<1>>{}) }
            -> std::same_as<StateVector<FiniteHilbertSpace<1>>>;
    };

    /// @brief Executor that schedules and runs a series of gate-like callables at compile-time (where possible).
    /// @tparam QBitCount  Number of qubits in the circuit.
    /// @tparam Gates      Variadic pack of gate-like callables accepted by the executor.
    template<natural_t QBitCount, QuantumGateLike... Gates>
    class QuantumCircuitExecutor
    {

        /// @brief Precompute 2^QBitCount for convenience
        static constexpr natural_t BasisStateCount = ConstexprMath::pow2(QBitCount);

        /// @brief The internal global state vector (amplitudes for 2^QBitCount basis states).
        StateVector<FiniteHilbertSpace<BasisStateCount>> m_stateVector;


        /// @brief Construct executor and immediately execute provided gates.
        /// @param gates  Variadic list of gate-like callables to apply in order.
        constexpr QuantumCircuitExecutor(const Gates& ... gates)
        {
            // Initialize to the |0...0> computational basis state
			m_stateVector = { complex_t::fromReal(1.0) };

            // Apply the provided gates in sequence
            executeCircuit(gates...);
        }

        friend class QuantumCircuit<QBitCount>;

    public:
        /// @brief Recursively apply gates: head then recurse on tail.
        /// @tparam Gate  First gate type.
        /// @tparam Rest  Remaining gate types.
        /// @param gate   Gate callable to apply.
        /// @param rest   Remaining gate callables.
        template<QuantumGateLike Gate, QuantumGateLike... Rest>
        constexpr void executeCircuit(const Gate& gate, const Rest&... rest)
        {
            // Apply the gate to the current state vector (gate returns a new vector)
            m_stateVector = gate(m_stateVector);

            // Recurse for the remaining gates
            executeCircuit(rest...);
        }

        /// @brief Base case for recursion: no gates left to apply.
        constexpr void executeCircuit() {}

        /// @brief Get the final state vector after executing all gates.
        constexpr const StateVector<FiniteHilbertSpace<BasisStateCount>>& getStateVector() const noexcept
        {
            return m_stateVector;
        }
    };

    /// @brief Facade class to create executors bound to a fixed qubit count.
    template<natural_t QBitCount>
    class QuantumCircuit
    {
    public:
        /// @brief Create an executor with the provided gate sequence.
        /// @tparam Gates  Gate-like callables to include in the circuit.
        /// @param gates   Instances of the gate-like callables (passed by reference-to-const).
        /// @return        A `QuantumCircuitExecutor` that has already executed the gates on an initialized state.
        template<QuantumGateLike... Gates>
        constexpr QuantumCircuitExecutor<QBitCount, Gates...> withGates(const Gates& ... gates) const
        {
            return QuantumCircuitExecutor<QBitCount, Gates...>(gates...);
        }
    };
}