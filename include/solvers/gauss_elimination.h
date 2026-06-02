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
        template<hilbert_space_t HilbertSpace, natural_t LevelCount>
        static constexpr StateVector<HilbertSpace> solve(const five_band_matrix_t<LevelCount>& FiveBandM, StateVector<HilbertSpace>& psi) noexcept
        {
            // 1. Initialize the full dense matrix  
			static constexpr natural_t SystemLevelCount = std::remove_cvref_t<decltype(FiveBandM)>::Dim;
            square_matrix_t<SystemLevelCount> M{};

            // Copy the banded structure into the dense NxN matrix
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                M[k][k] = FiveBandM[MAINDIAGONAL][k];

                // Superdiagonal (forward jump within the block)
                if (k + 1 < SystemLevelCount && (k % LevelCount) < LevelCount - 1)
                {
                    M[k][k + 1] = FiveBandM[SUPERDIAGONAL][k];
                }
                // Subdiagonal (backward jump within the block)
                if (k > 0 && (k % LevelCount) > 0)
                {
                    M[k][k - 1] = FiveBandM[SUBDIAGONAL][k];
                }
                // Upper far band (jump to the next block)
                if (k + LevelCount < SystemLevelCount)
                {
                    M[k][k + LevelCount] = FiveBandM[UPPER_FAR][k];
                }
                // Lower far band (jump to the previous block)
                if (k >= LevelCount)
                {
                    M[k][k - LevelCount] = FiveBandM[LOWER_FAR][k];
                }
            }

            // 2. Gaussian elimination with partial pivoting
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                // -- FIND PIVOT --
                // To ensure safe constexpr execution, we evaluate the squared magnitude (re^2 + im^2)
                natural_t pivot_row = k;
                real_t max_mag_sq = (M[k][k].re * M[k][k].re) + (M[k][k].im * M[k][k].im);

                for (natural_t i = k + 1; i < SystemLevelCount; ++i)
                {
                    real_t current_mag_sq = (M[i][k].re * M[i][k].re) + (M[i][k].im * M[i][k].im);

                    if (current_mag_sq > max_mag_sq)
                    {
                        max_mag_sq = current_mag_sq;
                        pivot_row = i;
                    }
                }

                // -- SWAP ROWS --
                // If a numerically larger element is found, swap the rows
                if (pivot_row != k)
                {
                    for (natural_t j = k; j < SystemLevelCount; ++j)
                    {
                        complex_t temp = M[k][j];
                        M[k][j] = M[pivot_row][j];
                        M[pivot_row][j] = temp;
                    }
                    // Swap the corresponding elements in the right-hand side vector (psi)
                    complex_t temp_psi = psi[k];
                    psi[k] = psi[pivot_row];
                    psi[pivot_row] = temp_psi;
                }

                const complex_t Pivot = M[k][k];

                // Safeguard against singular matrices (physically unlikely, but prevents division by zero)
                if (Pivot.re == 0.0 && Pivot.im == 0.0)
                {
                    continue;
                }

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