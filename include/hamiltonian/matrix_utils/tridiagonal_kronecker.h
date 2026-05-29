#pragma once
#include "core_types.h"
#include "solvers/matrix_access.h"

namespace KetCat
{
    /// @brief Computes the Kronecker (tensor) product of two tridiagonal matrices
    ///        and returns the resulting compact pentadiagonal matrix representation.
    ///
    /// @details
    /// Given two tridiagonal matrices A and B:
    ///
    ///     C = A ⊗ B
    ///
    /// the resulting matrix has a total flattened dimension of:
    ///
    ///     ResultDim = DimA * DimB
    ///
    /// The tensor product maps physical systems (such as a 2-qubit register) 
    /// into a joint Hilbert space. For independent 1D tridiagonal sub-operators 
    /// representing a 2D Laplace operator, the resulting Kronecker structure is 
    /// rigorously pentadiagonal. This function bypasses heavy dense matrix allocations 
    /// and populates the compact pentadiagonal_matrix_t framework directly.
    ///
    /// The flattened basis ordering follows:
    ///
    ///     |i,n⟩ -> i * DimB + n
    ///
    /// where:
    ///     i : state index in tridiagonal sub-system A
    ///     n : state index in tridiagonal sub-system B
    ///
    /// @tparam DimA Dimension of the first tridiagonal matrix.
    /// @tparam DimB Dimension of the second tridiagonal matrix.
    ///
    /// @param A First tridiagonal matrix operator.
    /// @param B Second tridiagonal matrix operator.
    ///
    /// @return Compact pentadiagonal tensor-product matrix representation.
    template<natural_t DimA, natural_t DimB>
    constexpr pentadiagonal_matrix_t<DimA* DimB>
        tensorProduct(
            const tridiagonal_matrix_t<DimA>& A,
            const tridiagonal_matrix_t<DimB>& B) noexcept
    {
        constexpr natural_t ResultDim = DimA * DimB;

        // Initialize the compact pentadiagonal structure with complex zeros
        pentadiagonal_matrix_t<ResultDim> Result{};

        // Setup internal types for structural access
        using AccessA = MatrixAccess<tridiagonal_matrix_t<DimA>>;
        using AccessB = MatrixAccess<tridiagonal_matrix_t<DimB>>;
        using AccessResult = MatrixAccess<pentadiagonal_matrix_t<ResultDim>>;

        // Multi-dimensional grid loop over subsystem tensor fields
        for (natural_t i = 0; i < DimA; ++i)
        {
            for (natural_t j = 0; j < DimA; ++j)
            {
                const complex_t ContribA = AccessA::get(A, i, j);

                // Skip unpopulated source elements to maximize loop efficiency
                if (ContribA == complex_t::zero())
                {
                    continue;
                }

                for (natural_t n = 0; n < DimB; ++n)
                {
                    for (natural_t m = 0; m < DimB; ++m)
                    {
                        const complex_t ContribB = AccessB::get(B, n, m);

                        // Skip unpopulated source elements
                        if (ContribB == complex_t::zero())
                        {
                            continue;
                        }

                        // Map multi-index structures into flattened global coordinates
                        const natural_t Row = i * DimB + n;
                        const natural_t Col = j * DimB + m;

                        // Transprently route the element via the unified matrix access trait
                        AccessResult::set(Result, Row, Col, ContribA * ContribB);
                    }
                }
            }
        }

        return Result;
    }
}
