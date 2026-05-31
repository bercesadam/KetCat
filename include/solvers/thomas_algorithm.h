#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    /// @brief Thomas algorithm backend.
    template<natural_t LevelCount>
    struct LinearSolver<LinearSolverBackend::ThomasTridiagonal, LevelCount>
    {
        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace>
            solve(tridiagonal_matrix_t<LevelCount> M,
                StateVector<HilbertSpace>& psi) noexcept
        {
            // Forward elimination
            for (natural_t i = 1; i < LevelCount; ++i)
            {
				// Compute the multiplier for the current row
                const complex_t Multiplier = M[SUBDIAGONAL][i] / M[MAINDIAGONAL][i - 1];

                M[MAINDIAGONAL][i] = M[MAINDIAGONAL][i] - Multiplier * M[SUPERDIAGONAL][i - 1];

                psi[i] = psi[i] - Multiplier * psi[i - 1];
            }

            // Back substitution
            StateVector<HilbertSpace> Result{};

            Result[LevelCount - 1] =  psi[LevelCount - 1] / M[MAINDIAGONAL][LevelCount - 1];

            for (natural_t i = LevelCount - 1; i-- > 0;)
            {
                Result[i] = (psi[i] - M[SUPERDIAGONAL][i] * Result[i + 1]) / M[MAINDIAGONAL][i];
            }

            return Result;
        }
    };
}
