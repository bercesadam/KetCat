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
    
    private:
        /// Embedded basis states |φ_i⟩ represented in the full grid space
        std::array<StateVector<FullHilbertSpace>, LevelCount> m_Basis;

        /// Energies of the basis states, computed as ⟨φ_i|H|φ_i⟩, where H is the system Hamiltonian.
        std::array<real_t, LevelCount> m_Energies;

    public:
        /// @brief Construct the reduced space from an arbitrary wavefunction generator.
        ///
        /// The generator must satisfy the `wavefunction_generator_t` concept.
        /// Each parameter tuple provides inputs to generate one basis state.
        ///
        /// @param generator Callable that produces a wavefunction from a parameter tuple
        /// @param params    Array of LevelCount parameter tuples for the generator
        template<typename Generator, typename ParamTuple>
            requires wavefunction_generator_t<Generator, ParamTuple>
        ReducedEnergySpace(Generator&& generator, const std::array<ParamTuple, LevelCount>& params)
        {
            for (natural_t i = 0; i < LevelCount; ++i)
            {
               auto w = std::apply(generator, params[i]);
               m_Basis[i] = w.m_Psi;
               m_Energies[i] = w.m_Energy;
            }

            orthonormalize();
        }


        /// @brief Project a full spatial state into the reduced subspace.
        ///
        /// Computes coefficients:
        ///
        ///      c_i = ⟨φ_i | ψ⟩
        ///
        /// Resulting reduced state:
        ///
        ///      |ψ⟩ → (c₀, ..., c_{K-1})
        constexpr StateVector<ReducedHilbertSpace>
            project(const StateVector<FullHilbertSpace>& psi) const noexcept
        {
            StateVector<ReducedHilbertSpace> Coeffs{};

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                Coeffs[i] = m_Basis[i].innerProduct(psi);
            }

            return Coeffs;
        }

        /// @brief Embed a reduced state back into the full grid space.
        ///
        /// Computes:
        ///
        ///      |ψ_red⟩ = Σ_i c_i |φ_i⟩
        constexpr StateVector<FullHilbertSpace>
            embed(const StateVector<ReducedHilbertSpace>& coeffs) const noexcept
        {
            StateVector<FullHilbertSpace> Psi{ complex_t::zero() };

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                for (natural_t k = 0; k < FullHilbertSpace::Dim; ++k)
                {
                    Psi[k] += coeffs[i] * m_Basis[i][k];
                }
            }

            return Psi;
        }

        /// @brief Get the energies of the reduced basis states.
        /// @return Array of energies corresponding to the basis states |φ_i⟩
        constexpr std::array<real_t, LevelCount> getEnergies() const noexcept
        {
            return m_Energies;
        }

    private:
        /// @brief Orthonormalize the basis using Gram–Schmidt.
        ///
        /// Ensures the basis satisfies:
        ///
        ///      ⟨φ_i | φ_j⟩ = δ_ij
        constexpr void orthonormalize() noexcept
        {
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                for (natural_t j = 0; j < i; ++j)
                {
                    auto Proj = m_Basis[j].innerProduct(m_Basis[i]);
                    for (natural_t k = 0; k < FullHilbertSpace::Dim; ++k)
                    {
                        m_Basis[i][k] = m_Basis[i][k] - Proj * m_Basis[j][k];
                    }
                }
                m_Basis[i].normalize();
            }
        }
    };

}