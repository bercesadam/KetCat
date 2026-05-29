#pragma once
#include "backend_traits.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    /// @brief Pentadiagonal Matrix Algorithm (PDMA) backend.
    ///
    /// @details
    /// This struct implements a specialized version of Gaussian elimination
    /// designed to solve linear systems with a pentadiagonal coefficient matrix
    ///
    ///   M · Result = psi
    ///
    /// in highly efficient O(N) execution time. This algorithm is an extension
    /// of the Thomas algorithm (which is restricted to tridiagonal systems).
    ///
    /// The pentadiagonal structure is a direct consequence of the 2D Laplace operator. When a 
    /// multi-dimensional system's Hamiltonian is constructed via the tensor product 
    /// (Kronecker product) of two independent 1D tridiagonal Hamiltonians, the 
    /// resulting joint operator exhibits a sparse, banded pentadiagonal form, which is used for
    /// solving two-atom systems (eg. CZ gates).
    ///
    /// The solver processes the pentadiagonal matrix M via a structured
    /// forward elimination pass to reduce it to an upper triangular form,
    /// followed by a back substitution pass to evaluate the solution vector.
    ///
    /// @tparam Dim Dimension of the linear system and the Hilbert space.
    template<natural_t Dim>
    struct LinearSolver<LinearSolverBackend::Pentadiagonal, Dim>
    {
        /// @brief  Solves the linear system M · Result = psi for a pentadiagonal matrix.
        /// @tparam HilbertSpace The underlying quantum Hilbert space structure.
        /// @param  M            The pentadiagonal coefficient matrix (passed by value to allow in-place modification).
        /// @param  psi          The right-hand side state vector (modified in-place).
        /// @return              The computed solution StateVector.
        ///
        /// @details
        /// The implementation assumes the matrix elements are packed into specific
        /// diagonal rows (SUBDIAGONAL2, SUBDIAGONAL, MAINDIAGONAL, SUPERDIAGONAL, SUPERDIAGONAL2).
        /// It performs sequential elimination and back substitution without pivoting,
        /// assuming the system is well-conditioned.
        template<hilbert_space_t HilbertSpace>
        static constexpr StateVector<HilbertSpace> solve(
            pentadiagonal_matrix_t<Dim> M,
            StateVector<HilbertSpace>& psi) noexcept
        {
            // Forward elimination
            for (natural_t i = 1; i < Dim; ++i)
            {
                // Eliminate the furthest lower diagonal (SUBDIAGONAL2) using row (i-2)
                if (i > 1)
                {
                    const complex_t Multiplier2 = M[SUBDIAGONAL2][i] / M[MAINDIAGONAL][i - 2];

                    M[SUBDIAGONAL][i] = M[SUBDIAGONAL][i] - Multiplier2 * M[SUPERDIAGONAL][i - 2];
                    M[MAINDIAGONAL][i] = M[MAINDIAGONAL][i] - Multiplier2 * M[SUPERDIAGONAL2][i - 2];

                    psi[i] = psi[i] - Multiplier2 * psi[i - 2];
                }

                // Eliminate the immediate lower diagonal (SUBDIAGONAL) using row (i-1)
                const complex_t Multiplier1 = M[SUBDIAGONAL][i] / M[MAINDIAGONAL][i - 1];

                M[MAINDIAGONAL][i] = M[MAINDIAGONAL][i] - Multiplier1 * M[SUPERDIAGONAL][i - 1];
                M[SUPERDIAGONAL][i] = M[SUPERDIAGONAL][i] - Multiplier1 * M[SUPERDIAGONAL2][i - 1];

                psi[i] = psi[i] - Multiplier1 * psi[i - 1];
            }

            // Back substitution
            StateVector<HilbertSpace> Result{};

            // Handle the boundary conditions at the bottom of the matrix
            Result[Dim - 1] = psi[Dim - 1] / M[MAINDIAGONAL][Dim - 1];

            if (Dim > 1)
            {
                Result[Dim - 2] = (psi[Dim - 2] - M[SUPERDIAGONAL][Dim - 2] * Result[Dim - 1])
                    / M[MAINDIAGONAL][Dim - 2];
            }

            // Back substitution loop for the remaining vector elements
            for (natural_t i = Dim - 2; i-- > 0; )
            {
                Result[i] = (psi[i] - M[SUPERDIAGONAL][i] * Result[i + 1]
                    - M[SUPERDIAGONAL2][i] * Result[i + 2])
                    / M[MAINDIAGONAL][i];
            }

            return Result;
        }
    };
}
