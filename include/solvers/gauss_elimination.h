#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
	/// @brief Dense Gaussian elimination backend - currently unused and superseded by the Pentadiagonal solvers,
    /// but retained for potential future use in non-banded systems.
    template<natural_t Dim>
    struct LinearSolver<LinearSolverBackend::GaussianElimination, Dim>
    {
        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace>
            solve(square_matrix_t<Dim> M,
                StateVector<HilbertSpace>& psi) noexcept
        {
            // Forward elimination
            for (natural_t k = 0; k < Dim; ++k)
            {
                const complex_t Pivot = M[k][k];

                for (natural_t j = k; j < Dim; ++j)
                {
                    M[k][j] = M[k][j] / Pivot;
                }

                psi[k] = psi[k] / Pivot;

                for (natural_t i = k + 1; i < Dim; ++i)
                {
                    const complex_t Factor = M[i][k];

                    for (natural_t j = k; j < Dim; ++j)
                    {
                        M[i][j] = M[i][j] - Factor * M[k][j];
                    }

                    psi[i] = psi[i] - Factor * psi[k];
                }
            }

            // Back substitution
            StateVector<HilbertSpace> Result{};

            for (natural_t i = Dim; i-- > 0;)
            {
                complex_t Sum = psi[i];

                for (natural_t j = i + 1; j < Dim; ++j)
                {
                    Sum = Sum - M[i][j] * Result[j];
                }

                Result[i] = Sum;
            }

            return Result;
        }
    };
}
