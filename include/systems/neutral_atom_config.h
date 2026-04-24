#pragma once
#include "core_types.h"

#include "hilbert_space/hilbert.h"

#include "atomic_physics_core/elements.h"
#include "atomic_physics_core/quantum_number.h"


namespace KetCat
{   
    /// @brief Compile-time configuration describing a neutral alkali-atom quantum computer's
    ///        fundamental properties within the KetCat framework
    ///
    /// @details
    /// This structure encodes all atom-specific, level-specific, and discretization
    /// parameters at compile time. It serves as a central configuration object for
    /// neutral-atom simulations and quantum-control models.
    ///
    /// The configuration specifies:
    ///  - the chemical element
    ///  - the spatial discretization of the wavefunction Hilbert space
    ///  - the logical and Rydberg levels used by higher-level abstractions
    ///  - the complete set of participating quantum numbers
    ///
    /// The design is fully static and type-driven, enabling strong compile-time
    /// guarantees and eliminating runtime configuration overhead.
    ///
    /// @tparam E
    ///     Chemical element of the atom (e.g. Cs, Rb, Na).
    /// @tparam Steps
    ///     Number of discretization points per spatial dimension for wavefunctions.
    /// @tparam Extent
    ///     Physical spatial extent of the wavefunction domain.
    /// @tparam Logical0
    ///     Index of the quantum level representing logical |0⟩.
    /// @tparam Logical1
    ///     Index of the quantum level representing logical |1⟩.
    /// @tparam _RydbergLevel
    ///     Index of the Rydberg (excited) level.
    /// @tparam _QuantumNumbers
    ///     Variadic list of quantum-number types defining the atomic level manifold.
    ///
    /// @note
    /// Logical and Rydberg levels are validated at compile time to ensure
    /// a consistent ordering within the supplied quantum-number list.
    ///
    template
        <
            Element E,
            natural_t Steps, real_t Extent,
            natural_t Logical0, natural_t Logical1, natural_t _RydbergLevel,
            quantum_number_t... _QuantumNumbers
        >
    struct NeutralAtomTypeConfig
    {
        template<DimensionTag SpatialDimensions>
        using HilbertSpaceStub = InfiniteHilbertSpace<SpatialDimensions, Steps, Extent>;

        static constexpr Element ChemicalElement = E;
        static constexpr natural_t WavefunctionDiscretizationSteps = Steps;
        static constexpr real_t WavefunctionPhysicalExtent = Extent;
        
        static constexpr natural_t Logical0Level = Logical0;
        static constexpr natural_t Logical1Level = Logical1;
        static constexpr natural_t RydbergLevel  = _RydbergLevel;

        static constexpr natural_t LevelCount = sizeof...(_QuantumNumbers);
        using QuantumNumbers = std::tuple<_QuantumNumbers...>;

        static_assert(LevelCount >= 3);
        static_assert(LevelCount - 3 >= Logical0);
        static_assert(LevelCount - 2 >= Logical1);
        static_assert(LevelCount - 1 >= RydbergLevel);
        static_assert(Logical0 < Logical1);
        static_assert(Logical1 < RydbergLevel);
    };
}

