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


    
    // @brief Concept for wavefunction generators constrained to a specific Hilbert space.
    /// @tparam Generator     Type of the wavefunction generator.
    /// @tparam ParamTuple    Tuple of parameters passed to the generator.
    /// @tparam HilbertSpace  The exact Hilbert space the generator must produce.
    template<
        typename Generator,
        typename ParamTuple,
        typename HilbertSpace
    >
    concept wavefunction_generator_for_t =
        hilbert_space_t<HilbertSpace> &&
        requires(Generator g, ParamTuple t)
    {
        { std::apply(g, t) } -> std::same_as<Wavefunction<HilbertSpace>>;
    };


	/// @brief Type alias for a fixed-size array of wavefunctions forming a basis set.
	/// @tparam HilbertSpace The type of the underlying Hilbert space for each wavefunction.
	/// @tparam NumStates    The number of wavefunctions in the basis set.
    template<hilbert_space_t HilbertSpace, natural_t NumStates>
    using basis_set_t = std::array<Wavefunction<HilbertSpace>, NumStates>;
}
