#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"


namespace KetCat
{
    /// @brief Computes the radial dipole transition matrix elements for a given basis set.
    ///
    /// @details
    ///   This function evaluates the transition matrix elements μₙₘ of the dipole operator 
    ///   acting on reduced radial wavefunctions u(r) = r·R(r). In the 1D radial representation, 
    ///   the integral simplifies to:
    ///
    ///   μₙₘ = ⟨ψₙ | r | ψₘ⟩ = ∫₀ᵟ uₙ*(r) · r · uₘ(r) dr
    ///
    ///   Since the dipole operator is Hermitian, only the upper triangle is explicitly computed, 
    ///   and the lower triangle is populated via complex conjugation: μₘₙ = μₙₘ*.
    ///
    /// @tparam HilbertSpace  A 1D spatial Hilbert space representing the radial grid [0, R_max].
    /// @tparam NumStates     The number of basis states to include in the matrix.
    ///
    /// @param basisStates    A collection of normalized radial states (e.g., 6s, 7p, 10g).
    /// @return               A complex-valued NumStates x NumStates dipole matrix.
    template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace, natural_t NumStates>
    constexpr matrix_t<NumStates> buildRadialDipoleMatrix(
        const basis_set_t<HilbertSpace, NumStates>& basisStates) noexcept
    {
        matrix_t<NumStates> DipoleMatrix{};

        constexpr natural_t Steps = HilbertSpace::Steps;
        constexpr real_t dx = HilbertSpace::dx();

        /// Iterate over state pairs (n, m) leveraging Hermiticity: μₙₘ = μₘₙ*
        for (natural_t n = 0; n < NumStates; ++n)
        {
            for (natural_t m = n; m < NumStates; ++m)
            {
                complex_t Integral = complex_t::zero();
                const auto& Bra = basisStates[n].m_Psi;
                const auto& Ket = basisStates[m].m_Psi;

                /// Numerical integration across the radial manifold:
                /// Integrand: uₙ*(r) · r · uₘ(r)
                /// The r² volume element is absorbed by the reduced radial normalization.
                for (natural_t k = 0; k < Steps; ++k)
                {
                    const real_t r = static_cast<real_t>(k) * dx;

                    /// Local transition density: uₙ*(rₖ) · rₖ · uₘ(rₖ)
                    const complex_t Term = Bra[k].conj() * r * Ket[k];
                    Integral = Integral + Term;
                }

                /// Apply the discrete integration measure Δr
                const complex_t MatrixElement = Integral * dx;

                /// Symmetry mapping for the off-diagonal elements
                DipoleMatrix[n][m] = MatrixElement;
                if (n != m)
                {
                    DipoleMatrix[m][n] = MatrixElement.conj();
                }
            }
        }

        return DipoleMatrix;
    }
}
