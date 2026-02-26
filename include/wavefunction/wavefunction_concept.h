#pragma once
#include "core_types.h"
#include "hilbert_space/hilbert.h"
    
namespace KetCat
{
    /// @brief Concept for wavefunction generators that produce StateVector objects in a specified Hilbert space.
    /// @tparam Generator  Type of the wavefunction generator (e.g. a struct with operator()).
    /// @tparam Args       Argument types required to generate the wavefunction (e.g. quantum numbers, grid parameters).
    template<typename Generator, typename ParamTuple>
    concept wavefunction_generator_t =
    requires(Generator g, ParamTuple t)
    {
        { std::apply(g, t) } -> std::same_as<
            StateVector<typename decltype(std::apply(g, t))::HilbertSpaceType>
        >;

        requires hilbert_space_t<typename decltype(std::apply(g, t))::HilbertSpaceType>;
    };
}
