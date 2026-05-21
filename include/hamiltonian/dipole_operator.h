#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"


namespace KetCat
{
    /// @brief Computes the full 3D dipole transition matrix elements leveraging analytical selection rules.
    ///
    /// @details
    ///   This function evaluates the transition matrix elements μₙₘ of the dipole operator (z = r·cosθ)
    ///   acting on the full 3D atomic states. Because the atomic potential is spherically symmetric,
    ///   the 3D wavefunction factors exactly into a radial part and an angular part:
    ///   Ψ(r, θ, φ) = R_nl(r) · Y_l^m(θ, φ).
    ///   
    ///   Consequently, the full 3D volume integral factors into a product of a radial integral 
    ///   and an angular integral:
    ///   ⟨ψₙ | z | ψₘ⟩ = ∫₀^R_max uₙ*(r)·r·uₘ(r) dr  ×  ∫∫ Y_lₙ^mₙ*(θ, φ) · cosθ · Y_lₘ^mₘ(θ, φ) sinθ dθ dφ
    ///
    ///   The angular integral evaluates analytically via the Wigner-Eckart theorem and spherical 
    ///   harmonic orthogonality, yielding strict selection rules: Δl = ±1 and Δm = 0.
    ///   Evaluating these selection rules via control flow is mathematically exact and avoids 
    ///   prohibitively expensive 3D numerical grid integration.
    ///
    /// @tparam HilbertSpace  A 1D spatial Hilbert space representing the radial grid [0, R_max].
    /// @tparam NumStates      The number of basis states to include in the matrix.
    ///
    /// @param basisStates    A collection of orthonormalized basis states
    /// @param quantumNumbers A array of pairs of quantum numbers (in l, m order) corresponding to each basis state.
    /// @return               A complex-valued NumStates x NumStates dipole matrix.
    template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace, natural_t NumStates>
    constexpr square_matrix_t<NumStates> calculateDipoleMatrix(
        const basis_set_t<HilbertSpace, NumStates>& basisStates,
        const std::array<std::pair<natural_t, natural_t>, NumStates>& quantumNumbers) noexcept
    {
        square_matrix_t<NumStates> DipoleMatrix{};

        constexpr natural_t Steps = HilbertSpace::Steps;
        constexpr real_t dx = HilbertSpace::dx();

        /// Iterate over state pairs (n, m) leveraging Hermiticity: μₙₘ = μₘₙ*
        for (natural_t n = 0; n < NumStates; ++n)
        {
            const auto& BraState = basisStates[n];
			const natural_t BraL = quantumNumbers[n].first;
			const natural_t BraM = quantumNumbers[n].second;

            for (natural_t m = n; m < NumStates; ++m)
            {
                const auto& KetState = basisStates[m];
				const natural_t KetL = quantumNumbers[m].first;
				const natural_t KetM = quantumNumbers[m].second;

                /// Enforce angular selection rules derived from the spherical harmonics integral:
                /// 1. Parity and angular momentum conservation require Δl = l_bra - l_ket = ±1.
                /// 2. For a z-polarized field (operator ∝ cosθ), conservation of the projection requires Δm = m_bra - m_ket = 0.
                /// If these conditions are violated, the spatial angular integral is strictly zero.
				const int DeltaL = BraL - KetL;
                const int DeltaM = BraM - KetM;

                if (ConstexprMath::abs(DeltaL) != 1 || DeltaM != 0)
                {
                    /// Elements default to complex_t::zero() via matrix initialization.
                    /// Short-circuiting here bypasses the O(Steps) radial integration loop.
                    continue;
                }

                complex_t Integral = complex_t::zero();
                const auto& Bra = BraState.m_Psi;
                const auto& Ket = KetState.m_Psi;

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

                /// Analytical evaluation of the angular overlap integral under z-polarization:
                /// AngularFactor = ⟨l_bra, m_bra | cosθ | l_ket, m_ket⟩
                /// Computed via standard 3-j symbols or algebraic reduction.
                const real_t lMax = static_cast<real_t>(std::max(BraL, KetL));
                const real_t lBraSq = static_cast<real_t>(BraL * BraL);
                const real_t mBraSq = static_cast<real_t>(BraM * BraM);

                const real_t AngularFactor = std::sqrt((lMax * lMax - mBraSq) / (4.0 * lMax * lMax - 1.0));

                /// Combine the numerical radial matrix element with the analytical angular component,
                /// scaled by the discrete integration measure Δr
                const complex_t MatrixElement = Integral * dx * AngularFactor;

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