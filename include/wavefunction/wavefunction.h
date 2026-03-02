#pragma once
#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "hilbert_space/state_vector.h"

namespace KetCat
{
    /// @brief Wavefunction struct to hold a state vector and its associated energy.
    /// @tparam _HilbertSpace The type of the underlying Hilbert space, which determines the structure of the state vector.
    template <spatial_hilbert_space_t _HilbertSpace>
    struct Wavefunction
    {
        /// @brief Expose the underlying Hilbert space type
        using HilbertSpace = _HilbertSpace;

        /// @brief State vector representing the wavefunction in the specified Hilbert space.
        StateVector<HilbertSpace> m_Psi;

        /// @brief  Hartree energy of the wavefunction, computed as ⟨ψ|H|ψ⟩.
        real_t m_Energy;
    };

    /// @brief Concept for wavefunction generators that produce `Wavefunction` objects in a specified Hilbert space.
    /// @tparam Generator  Type of the wavefunction generator (e.g. a struct with operator()).
    /// @tparam ParamTuple Argument types required to generate the wavefunction (e.g. quantum numbers, grid parameters).
    template<typename Generator, typename ParamTuple>
    concept wavefunction_generator_t =
        requires(Generator g, ParamTuple t)
    {
        // The generator must return a Wavefunction object (not just the raw state
        // vector) so that reduced-space constructors can access both the vector and
        // its associated energy.
        { std::apply(g, t) } -> std::same_as<
              Wavefunction<typename decltype(std::apply(g, t))::HilbertSpace>
              >;

        // Ensure the generated object refers to a valid Hilbert space type.
        requires hilbert_space_t<typename decltype(std::apply(g, t))::HilbertSpace>;
    };

}
