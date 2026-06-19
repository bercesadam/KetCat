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
    template<hilbert_space_t HilbertSpace>
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
        static StateVector<HilbertSpace, QuantumPicture::Dirac> toDiracPicture(
            StateVector<HilbertSpace, QuantumPicture::Schrodinger>& psi,
            const std::array<real_t, HilbertSpace::Dim>& diagEnergies) noexcept
        {
            StateVector<HilbertSpace, QuantumPicture::Dirac> Result{};
            Result.m_TimeStamp = psi.m_TimeStamp;

            for (natural_t i = 0; i < HilbertSpace::Dim; ++i)
            {
                const real_t PhaseArg = diagEnergies[i] * psi.m_TimeStamp;

                const complex_t Phase =
                {
                    std::cos(PhaseArg),
                    std::sin(PhaseArg)   // +i
                };

                Result[i] = psi[i] * Phase;
            }

            return Result;
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
        static StateVector<HilbertSpace, QuantumPicture::Schrodinger> toSchrodingerPicture(
            StateVector<HilbertSpace, QuantumPicture::Dirac>& psi,
            const std::array<real_t, HilbertSpace::Dim>& diagEnergies) noexcept
        {
            StateVector<HilbertSpace, QuantumPicture::Schrodinger> Result{};
            Result.m_TimeStamp = psi.m_TimeStamp;

            for (natural_t i = 0; i < HilbertSpace::Dim; ++i)
            {
                const real_t PhaseArg = diagEnergies[i] * psi.m_TimeStamp;

                const complex_t Phase =
                {
                    std::cos(PhaseArg),
                    -std::sin(PhaseArg)   // -i
                };

                Result[i] = psi[i] * Phase;
            }

            return Result;
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
    template<natural_t LevelCount, natural_t QubitCount>
    class InteractionPictureHamiltonian
    {
        /// @brief Local two-qubit subspace dimension ($LevelCount \times LevelCount$, e.g., 36).
        static constexpr natural_t Dim = ConstexprMath::pow(LevelCount, QubitCount);

    public:
        /**
            * @brief Transforms the Schrödinger-picture (RWA) 5-band matrix into the Interaction Picture.
            *
            * @param H
            * The current RWA Hamiltonian matrix containing Detuning, AC Stark shifts, and vdW blockade.
            * @param localRwaEnergies
            * The pre-calculated pure RWA detuning baseline array for the 2-qubit subspace (size: Dim).
            * @param t
            * The current absolute global simulation time.
            *
            * @return matrix_t
            * The effective Interaction Picture Hamiltonian matrix ready for the numerical solver.
            */
        template<BandedMatrix MatrixType>
        static MatrixType transform(const MatrixType& H,
            const eigenenergies_t<Dim>& localRwaEnergies,
            const real_t t) noexcept
        {
            MatrixType HI{};

            // --- 1. MAIN DIAGONAL: Strip the pure RWA detuning baseline ---
            // Since H[MAINDIAGONAL][i] contains (Detuning + StarkShift + vdW) and localRwaEnergies[i]
            // contains only the pure Detuning, the subtraction leaves the dynamic Stark shifts and 
            // static vdW blockade energies completely intact in the main diagonal.
            for (natural_t i = 0; i < Dim; ++i)
            {
                HI[MAINDIAGONAL][i] = H[MAINDIAGONAL][i] - complex_t::fromReal(localRwaEnergies[i]);
            }

            // --- 2. SUPER DIAGONAL (i → i+1) ---
            // Coherent couplings are rotated using the pure RWA detuning energy differences.
            for (natural_t i = 0; i < Dim - 1; ++i)
            {
                auto Hij = H[SUPERDIAGONAL][i];
                real_t dE = localRwaEnergies[i + 1] - localRwaEnergies[i];
                HI[SUPERDIAGONAL][i] = Hij * exp_i(dE * t);
            }

            // --- 3. SUB DIAGONAL (i → i-1) ---
            for (natural_t i = 1; i < Dim; ++i)
            {
                auto Hij = H[SUBDIAGONAL][i];
                real_t dE = localRwaEnergies[i - 1] - localRwaEnergies[i];
                HI[SUBDIAGONAL][i] = Hij * exp_i(dE * t);
            }

            if constexpr (FiveBandMatrix<MatrixType>)
            {
                // --- 4. FAR DIAGONALS (Cross-talk / Inter-qubit transitions) ---
                for (natural_t i = 0; i < Dim; ++i)
                {
                    if (i + LevelCount < Dim)
                    {
                        auto Hij = H[UPPER_FAR][i];
                        real_t dE = localRwaEnergies[i + LevelCount] - localRwaEnergies[i];
                        HI[UPPER_FAR][i] = Hij * exp_i(dE * t);
                    }

                    if (i >= LevelCount)
                    {
                        auto Hij = H[LOWER_FAR][i];
                        real_t dE = localRwaEnergies[i - LevelCount] - localRwaEnergies[i];
                        HI[LOWER_FAR][i] = Hij * exp_i(dE * t);
                    }
                }
            }

            return HI;
        }

    private:
        /**
            * @brief Evaluates the complex phase factor $e^{i x}$ using compile-time compatible math.
            */
        static complex_t exp_i(real_t x) noexcept
        {
            return {
                std::cos(x),
                std::sin(x)
            };
        }
    };
}
