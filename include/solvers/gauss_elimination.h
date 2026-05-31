#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    template<natural_t LevelCount>
    struct LinearSolver<LinearSolverBackend::FiveBandGaussianElimination, LevelCount>
    {
        static constexpr natural_t SystemLevelCount = LevelCount * LevelCount;

        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace>
            solve(const five_band_matrix_t<LevelCount>& FiveBandM,
                  StateVector<HilbertSpace>& psi) noexcept
        {
            square_matrix_t<SystemLevelCount> M{};

            for (natural_t i = 0; i < SystemLevelCount; ++i) {
                M[i][i] = FiveBandM[MAINDIAGONAL][i];
                natural_t n = i % LevelCount;
                if (i + 1 < SystemLevelCount && n < LevelCount - 1) {
                    M[i][i + 1] = FiveBandM[SUPERDIAGONAL][i];
                }
                if (i > 0 && n > 0) {
                    M[i][i - 1] = FiveBandM[SUBDIAGONAL][i];
                }
                if (i + LevelCount < SystemLevelCount) {
                    M[i][i + LevelCount] = FiveBandM[UPPER_FAR][i];
                }
                if (i >= LevelCount) {
                    M[i][i - LevelCount] = FiveBandM[LOWER_FAR][i];
                }
            }

            // Forward elimination
            for (natural_t k = 0; k < SystemLevelCount; ++k)
            {
                const complex_t Pivot = M[k][k];

                for (natural_t j = k; j < SystemLevelCount; ++j)
                {
                    M[k][j] = M[k][j] / Pivot;
                }

                psi[k] = psi[k] / Pivot;

                for (natural_t i = k + 1; i < SystemLevelCount; ++i)
                {
                    const complex_t Factor = M[i][k];

                    for (natural_t j = k; j < SystemLevelCount; ++j)
                    {
                        M[i][j] = M[i][j] - Factor * M[k][j];
                    }

                    psi[i] = psi[i] - Factor * psi[k];
                }
            }

            // Back substitution
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
