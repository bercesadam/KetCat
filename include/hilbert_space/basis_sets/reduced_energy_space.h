#pragma once
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"

namespace KetCat
{

    /// @brief Reduced Hilbert space constructed from an arbitrary set of wavefunction seeds.
    ///
    /// @tparam OriginalHilbertSpace The original full Hilbert space (e.g. an infinite-dimensional spatial grid).
    /// @tparam LevelCount Number of retained basis states in the reduced space.
    ///
    /// ------------------------------------------------------------------------
    /// Mathematical construction
    /// ------------------------------------------------------------------------
    ///
    /// Full spatial Hilbert space:
    ///
    ///      𝓗_grid  ≅  ℂ^(GridDim²)
    ///
    /// Seed wavefunctions:
    ///
    ///      ψ_i(x,y)  (provided by a user-defined generator)
    ///
    /// We select LevelCount states:
    ///
    ///      { |φ₀⟩, |φ₁⟩, ..., |φ_{K-1}⟩ }
    ///
    /// and define the reduced subspace:
    ///
    ///      𝓗_red = span{ |φ_i⟩ }
    ///
    /// which is isomorphic to ℂ^K.
    ///
    /// This provides a physically meaningful qudit or truncated subspace
    /// suitable for simulation.
    ///
    /// ------------------------------------------------------------------------
    /// Physical meaning
    /// ------------------------------------------------------------------------
    ///
    /// Each basis vector corresponds to a physical mode of the system
    /// as defined by the generator.
    ///
    /// Logical states can later be defined as selected combinations
    /// of these physical basis states.
    ///
    /// This object bridges:
    ///
    ///      InfiniteHilbertSpace ↔ FiniteHilbertSpace
    ///
    /// allowing full state-vector simulation in a reduced, physically
    /// meaningful subspace.
    template<spatial_hilbert_space_t InitialHilbertSpace, natural_t LevelCount>
    class ReducedEnergySpace
    {
    public:
        using FullHilbertSpace = InitialHilbertSpace;
        using ReducedHilbertSpace = FiniteHilbertSpace<LevelCount>;

        /// @brief Project a full spatial state into the reduced subspace.
        ///
        /// Computes coefficients:
        ///
        ///      c_i = ⟨φ_i | ψ⟩
        ///
        /// Resulting reduced state:
        ///
        ///      |ψ⟩ → (c₀, ..., c_{K-1})
		template <spatial_hilbert_space_t HilbertSpace>
        constexpr StateVector<ReducedHilbertSpace>
            project(const basis_set_t<HilbertSpace, LevelCount>& basisSet,
                    const StateVector<HilbertSpace>& psi) const noexcept
        {
            StateVector<ReducedHilbertSpace> Coeffs{};

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                Coeffs[i] = basisSet[i].m_Psi.innerProduct(psi);
            }

            return Coeffs;
        }

        /// @brief Embed a reduced state back into the full grid space.
        ///
        /// Computes:
        ///
        ///      |ψ_red⟩ = Σ_i c_i |φ_i⟩
        template <spatial_hilbert_space_t HilbertSpace>
        constexpr StateVector<HilbertSpace>
            embed(const basis_set_t<HilbertSpace, LevelCount>& basisSet,
                  const StateVector<ReducedHilbertSpace>& coeffs) const noexcept
        {
            StateVector<HilbertSpace> Psi{ complex_t::zero() };
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                for (natural_t k = 0; k < HilbertSpace::Dim; ++k)
                {
                    Psi[k] += coeffs[i] * basisSet[i].m_Psi[k];
                }
            }

            return Psi;
        }
    };

}