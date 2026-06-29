#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    /// @brief Direct linear system solver utilizing full dense representation populated from 2D tensor-product 5-band matrices.
    /// @details Performs standard Gaussian elimination by reconstructing the full matrix layout 
    /// originating from a 2D Hamiltonian. It maps the 5-band structure into a temporary dense 
    /// N × N representation to safely perform forward elimination and backward substitution.
    /// The global dimension is determined by N = LevelCount².
    /// @tparam LevelCount The number of single-atom eigenstates along a single dimension.
    template<natural_t LevelCount>
    struct LinearSolver<LinearSolverBackend::FiveBandGaussianElimination, LevelCount>
    {
        /// @brief Solves the linear system M · x = ψ using dense Gaussian elimination.
        /// @tparam HilbertSpace Type configuration of the underlying atomic state space.
        /// @tparam LevelCount The number of single-atom eigenstates.
        /// @param FiveBandM Compact 5-band representation containing the Hamiltonian operators.
        /// @param psi Input right-hand side state vector |ψ⟩ to be mutated during elimination.
        /// @return The computed state vector |x⟩ solving the linear system.
        template<hilbert_space_t HilbertSpace, natural_t LevelCount>
        static constexpr StateVector<HilbertSpace> solve(const five_band_matrix_t<LevelCount>& FiveBandM, StateVector<HilbertSpace>& psi) noexcept
        {
            // 1. Initialize the full dense matrix  
            static constexpr natural_t SystemLevelCount = std::remove_cvref_t<decltype(FiveBandM)>::Dim;
            square_matrix_t<SystemLevelCount> M{};

            // Copy the banded structure into the dense NxN matrix safely
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                M[k][k] = FiveBandM[MAINDIAGONAL][k];

                // SUPERDIAGONAL: k-th row, k+1-th column
                // Only copy if we do not cross the boundary of the block size defined by LevelCount
                if (k + 1 < SystemLevelCount && (k % LevelCount) < LevelCount - 1)
                {
                    M[k][k + 1] = FiveBandM[SUPERDIAGONAL][k];
                }

                // SUBDIAGONAL: k-th row, k-1-th column
                // Only copy if we are not at the first element inside the internal block boundary
                if (k > 0 && (k % LevelCount) > 0)
                {
                    M[k][k - 1] = FiveBandM[SUBDIAGONAL][k];
                }

                // UPPER_FAR: k-th row, k + LevelCount-th column
                // This is not constrained by the inner block index as it represents cross-block transitions
                if (k + LevelCount < SystemLevelCount)
                {
                    M[k][k + LevelCount] = FiveBandM[UPPER_FAR][k];
                }

                // LOWER_FAR: k-th row, k - LevelCount-th column
                // This is not constrained by the inner block index as it represents cross-block transitions
                if (k >= LevelCount)
                {
                    M[k][k - LevelCount] = FiveBandM[LOWER_FAR][k];
                }
            }

            // 2. Gaussian elimination
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                const complex_t Pivot = M[k][k];

                // -- NORMALIZE PIVOT ROW --
                for (natural_t j = k; j < SystemLevelCount; ++j)
                {
                    M[k][j] = M[k][j] / Pivot;
                }
                psi[k] = psi[k] / Pivot;

                // -- ELIMINATE ROWS BELOW --
                for (natural_t i = k + 1; i < SystemLevelCount; ++i)
                {
                    const complex_t Factor = M[i][k];

                    // Skip calculations if the factor is already exactly zero
                    if (Factor.re == 0.0 && Factor.im == 0.0)
                    {
                        continue;
                    }

                    for (natural_t j = k; j < SystemLevelCount; ++j)
                    {
                        M[i][j] = M[i][j] - Factor * M[k][j];
                    }
                    psi[i] = psi[i] - Factor * psi[k];
                }
            }

            // 3. Back substitution
            // Accurately determine the solution vector from the upper triangular matrix
            StateVector<HilbertSpace> Result{};
            for (natural_t i = SystemLevelCount; i-- > 0;)
            {
                complex_t Sum = psi[i];
                for (natural_t j = i + 1; j < SystemLevelCount; ++j)
                {
                    Sum = Sum - M[i][j] * Result[j];
                }
                Result[i] = Sum;
            }

            return Result;
        }
    };

}