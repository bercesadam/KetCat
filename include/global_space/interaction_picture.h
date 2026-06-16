#pragma once
#include "core_types.h"

namespace KetCat
{
    /// @brief State transformation between Schrödinger and Interaction (Dirac) pictures
    ///
    /// @details
    /// Handles transformations:
    ///
    ///   |ψ_I(t)⟩ = exp(+i H0 t) |ψ(t)⟩      (to Dirac / interaction)
    ///   |ψ(t)⟩   = exp(-i H0 t) |ψ_I(t)⟩    (to Schrödinger)
    ///
    /// Assumes:
    ///   • H0 is diagonal
    ///   • diagEnergies[i] = ⟨i|H0|i⟩
    ///
    /// Works on the full global Hilbert space/state vector
    ///
    template<natural_t HilbertDim>
    class InteractionPictureStateTransformer
    {
    public:

        /// @brief Transform Schrödinger → Interaction (Dirac) picture
        ///
        /// @param psi             Global state vector (modified in-place)
        /// @param diagEnergies    Diagonal energies E_i of H0
        /// @param t               Absolute time
        ///
        /// @details
        /// Applies:
        ///     ψ_I[i] = ψ[i] * exp(+i E_i t)
        ///
        template<typename StateVectorType>
        static void toDiracPicture(
            StateVectorType& psi,
            const std::array<real_t, HilbertDim>& diagEnergies,
            real_t t) noexcept
        {
            for (natural_t i = 0; i < HilbertDim; ++i)
            {
                const real_t PhaseArg = diagEnergies[i] * t;

                const complex_t Phase =
                {
                    ConstexprMath::cos(PhaseArg),
                    ConstexprMath::sin(PhaseArg)   // +i
                };

                psi[i] = psi[i] * Phase;
            }
        }

        /// @brief Transform Interaction (Dirac) → Schrödinger picture
        ///
        /// @param psi             Global state vector (modified in-place)
        /// @param diagEnergies    Diagonal energies E_i of H0
        /// @param t               Absolute time
        ///
        /// @details
        /// Applies:
        ///     ψ[i] = ψ_I[i] * exp(-i E_i t)
        ///
        template<typename StateVectorType>
        static void toSchrodingerPicture(
            StateVectorType& psi,
            const std::array<real_t, HilbertDim>& diagEnergies,
            real_t t) noexcept
        {
            for (natural_t i = 0; i < HilbertDim; ++i)
            {
                const real_t PhaseArg = diagEnergies[i] * t;

                const complex_t Phase =
                {
                    ConstexprMath::cos(PhaseArg),
                   -ConstexprMath::sin(PhaseArg)   // -i
                };

                psi[i] = psi[i] * Phase;
            }
        }
    };

    /// @brief Interaction-picture Hamiltonian transformer for banded 2-atom matrix
    ///
    /// @details
    /// Transforms:
    ///
    ///     H → H_I(t) = exp(+i H0 t) H exp(-i H0 t) - H0
    ///
    /// Assumes:
    ///   • H0 = diagonal part of H
    ///   • Only off-diagonals survive
    ///   • Uses banded structure (no dense ops!)
    ///
    template<natural_t LevelCount>
    class InteractionPictureHamiltonian
    {
        static constexpr natural_t Dim = LevelCount * LevelCount;

    public:
        using matrix_t = five_band_matrix_t<LevelCount>;

        static matrix_t transform(
            const matrix_t& H,
            real_t t) noexcept
        {
            matrix_t HI{};

            const auto& MainDiagonal = H[MAINDIAGONAL];

            // --- DIAGONAL → ZERO ---
            for (natural_t i = 0; i < Dim; ++i)
            {
                HI[MAINDIAGONAL][i] = complex_t::zero();
            }

            // --- SUPER DIAGONAL (i → i+1)
            for (natural_t i = 0; i < Dim - 1; ++i)
            {
                auto Hij = H[SUPERDIAGONAL][i];

                real_t dE = MainDiagonal[i].re - MainDiagonal[i + 1].re;
                HI[SUPERDIAGONAL][i] = Hij * exp_i(dE * t);
            }

            // --- SUB DIAGONAL (i → i-1)
            for (natural_t i = 1; i < Dim; ++i)
            {
                auto Hij = H[SUBDIAGONAL][i];

                real_t dE = MainDiagonal[i].re - MainDiagonal[i - 1].re;
                HI[SUBDIAGONAL][i] = Hij * exp_i(dE * t);
            }

            // --- FAR DIAGONALS
            for (natural_t i = 0; i < Dim; ++i)
            {
                if (i + LevelCount < Dim)
                {
                    auto Hij = H[UPPER_FAR][i];

                    real_t dE = MainDiagonal[i].re - MainDiagonal[i + LevelCount].re;
                    HI[UPPER_FAR][i] = Hij * exp_i(dE * t);
                }

                if (i >= LevelCount)
                {
                    auto Hij = H[LOWER_FAR][i];

                    real_t dE = MainDiagonal[i].re - MainDiagonal[i - LevelCount].re;
                    HI[LOWER_FAR][i] = Hij * exp_i(dE * t);
                }
            }

            return HI;
        }

    private:
        static complex_t exp_i(real_t x) noexcept
        {
            return {
                ConstexprMath::cos(x),
                ConstexprMath::sin(x)
            };
        }
    };
}
