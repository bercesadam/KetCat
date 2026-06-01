#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    /// @brief Direct linear system solver optimized for 2D tensor-product 5-band matrices.
    /// @details Performs Gaussian elimination by exploiting the sparse structure originating 
    /// from a 2D Hamiltonian. It avoids full dense matrix pre-allocation by streaming 
    /// coefficients from the 5-band representation on-the-fly, reducing the forward 
    /// elimination complexity from O(N³) to O(N · LevelCount²).
    /// @tparam LevelCount The number of spatial discretization levels along a single dimension.
    template<natural_t LevelCount>
    struct LinearSolver<LinearSolverBackend::FiveBandGaussianElimination, LevelCount>
    {
        static constexpr natural_t SystemLevelCount = LevelCount * LevelCount;

        /// @brief Solves the linear system M · psi_next = psi using optimized banded Gaussian elimination.
        /// @param FiveBandM Compact storage of the 5-band input matrix operator.
        /// @param psi Right-hand side state vector (mutated in-place during elimination).
        /// @return The solved state vector at the next time step.
        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace>
            solve(const five_band_matrix_t<LevelCount>& FiveBandM,
                StateVector<HilbertSpace>& psi) noexcept
        {
            // Dense upper-triangular target matrix for back substitution (zero-initialized)
            square_matrix_t<SystemLevelCount> M{};

            // Forward elimination with localized banded window execution
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                // Stream current row operator elements from 5-band representation
                M[k][k] = FiveBandM[MAINDIAGONAL][k];

                const natural_t n_k = k % LevelCount;
                if (k + 1 < SystemLevelCount && n_k < LevelCount - 1)
                {
                    M[k][k + 1] = FiveBandM[SUPERDIAGONAL][k];
                }
                if (k + LevelCount < SystemLevelCount)
                {
                    M[k][k + LevelCount] = FiveBandM[UPPER_FAR][k];
                }

                // Normalize the pivot row within the upper band boundary
                const complex_t Pivot = M[k][k];
                const natural_t max_col = (k + LevelCount < SystemLevelCount) ? (k + LevelCount + 1) : SystemLevelCount;

                for (natural_t j = k; j < max_col; ++j)
                {
                    M[k][j] = M[k][j] / Pivot;
                }
                psi[k] = psi[k] / Pivot;

                // Eliminate sub-diagonal entries within the lower band boundary
                const natural_t max_row = (k + LevelCount < SystemLevelCount) ? (k + LevelCount + 1) : SystemLevelCount;

                for (natural_t i = k + 1; i < max_row; ++i)
                {
                    // Fetch lower band element if not yet modified by fill-in processing
                    if (M[i][k] == complex_t{ 0.0, 0.0 })
                    {
                        if (i == k + 1 && (k % LevelCount) < LevelCount - 1) {
                            M[i][k] = FiveBandM[SUBDIAGONAL][i];
                        }
                        else if (i == k + LevelCount) {
                            M[i][k] = FiveBandM[LOWER_FAR][i];
                        }
                    }

                    const complex_t Factor = M[i][k];
                    if (Factor == complex_t{ 0.0, 0.0 }) continue;

                    // Compute row subtraction limited to the current active band window
                    for (natural_t j = k; j < max_col; ++j)
                    {
                        M[i][j] = M[i][j] - Factor * M[k][j];
                    }
                    psi[i] = psi[i] - Factor * psi[k];
                }
            }

            // Back substitution phase utilizing the populated upper-triangular matrix
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