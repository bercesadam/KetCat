#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    /// @brief Direct O(N) linear system solver utilizing the Thomas algorithm.
    /// @details Implements an optimized tridiagonal matrix solver (simplified Gaussian 
    /// elimination) designed for 1D spatial discretization grids. 
    /// @tparam LevelCount The number of spatial discretization levels (system size).
    template<natural_t LevelCount>
    struct LinearSolver<LinearSolverBackend::ThomasTridiagonal, LevelCount>
    {
        /// @brief Solves the linear system M · psi_next = psi using the Thomas algorithm.
        /// @param M Compact tridiagonal input matrix operator (passed by value to allow in-place modification).
        /// @param psi Right-hand side state vector (mutated in-place during elimination).
        /// @return The solved state vector at the next time step.
        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace>
            solve(tridiagonal_matrix_t<LevelCount> M,
                StateVector<HilbertSpace>& psi) noexcept
        {
            // Forward elimination: eliminate the subdiagonal entries
            for (natural_t i = 1; i < LevelCount; ++i)
            {
                const complex_t Multiplier = M[SUBDIAGONAL][i] / M[MAINDIAGONAL][i - 1];

                M[MAINDIAGONAL][i] = M[MAINDIAGONAL][i] - Multiplier * M[SUPERDIAGONAL][i - 1];
                psi[i] = psi[i] - Multiplier * psi[i - 1];
            }

            // Back substitution phase
            StateVector<HilbertSpace> Result{};

            // Boundary condition for the last grid point
            Result[LevelCount - 1] = psi[LevelCount - 1] / M[MAINDIAGONAL][LevelCount - 1];

            // Backward recurrence loop for the remaining states
            for (natural_t i = LevelCount - 1; i-- > 0;)
            {
                Result[i] = (psi[i] - M[SUPERDIAGONAL][i] * Result[i + 1]) / M[MAINDIAGONAL][i];
            }

            return Result;
        }
    };
}