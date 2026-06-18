#pragma once
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"
#include "matrix_utils/matrix.h"


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
    class LocalSpaceHelper
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
       /// 
       /// @warning This is a purely mathematical embedding that does not
       /// preserve phase relations under entanglement, function kept for legacy
       /// and debug purposes. For visualization of entangled states, use the
       /// density-matrix version below.
        template <spatial_hilbert_space_t HilbertSpace>
        constexpr StateVector<HilbertSpace, QuantumPicture::Schrodinger>
            embed(const basis_set_t<HilbertSpace, LevelCount>& basisSet,
                const StateVector<ReducedHilbertSpace, QuantumPicture::Schrodinger>& coeffs) const noexcept
        {
            StateVector<HilbertSpace, QuantumPicture::Schrodinger> Psi{ complex_t::zero() };

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                for (natural_t k = 0; k < HilbertSpace::Dim; ++k)
                {
                    Psi[k] += coeffs[i] * basisSet[i].m_Psi[k];
                }
            }

            return Psi;
        }

       /// @brief Extract a representative, phase-coherent state vector within the reduced subspace from a density matrix.
       /// @details Extracts the dominant coherent state vector components directly from the density matrix.
       ///          This avoids full spatial grid embedding, returning the complex amplitudes inside the 
       ///          FiniteHilbertSpace<LevelCount> while preserving essential phase relations under entanglement.
       /// @param rho The reduced density matrix of the qubit.
       /// @return    A normalized state vector in the reduced finite Hilbert space.
        constexpr StateVector<ReducedHilbertSpace, QuantumPicture::Schrodinger>
            extractCoherentState(const Matrix<LevelCount>& rho) const noexcept
        {
            // 1. Find the pivot column with the largest diagonal population.
            natural_t Pivot = 0;
            real_t MaxPop = rho.at(0, 0).re;
            for (natural_t i = 1; i < LevelCount; ++i)
            {
                if (rho.at(i, i).re > MaxPop)
                {
                    MaxPop = rho.at(i, i).re;
                    Pivot = i;
                }
            }

            // 2. Extract the complex coefficients from the dominant pivot column.
            StateVector<ReducedHilbertSpace, QuantumPicture::Schrodinger> Coeffs{};
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                Coeffs[i] = rho.at(i, Pivot);
            }

            // 3. Normalize the extracted coefficients to ensure unit norm.
            real_t Norm2 = real_t(0);
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                Norm2 += Coeffs[i].normSquared();
            }

            if (Norm2 > real_t(0))
            {
                const real_t InvNorm = real_t(1) / ConstexprMath::sqrt(Norm2);
                for (natural_t i = 0; i < LevelCount; ++i)
                {
                    Coeffs[i] = Coeffs[i] * InvNorm;
                }
            }

            return Coeffs;
        }

       
    };
}
