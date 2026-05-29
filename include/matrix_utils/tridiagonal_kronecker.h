#pragma once
#include "core_types.h"
#include "operation_space/utils/matrix.h"

namespace KetCat
{
    /// @brief Computes the Kronecker (tensor) product of two tridigonal matrices
    ///        and returns the resulting dense matrix representation.
    ///
    /// @details
    /// Given two tridigonal matrices A and B:
    ///
    ///     C = A ⊗ B
    ///
    /// the resulting matrix has dimension:
    ///
    ///     (DimA * DimB) × (DimA * DimB)
    ///
    /// The tensor product is expanded explicitly into a dense matrix_t,
    /// since the resulting structure is generally no longer tridigonal.
    ///
    /// The flattened basis ordering is:
    ///
    ///     |i,j⟩ -> i * DimB + j
    ///
    /// where:
    ///     i : state index in A
    ///     j : state index in B
    ///
    /// @tparam DimA Dimension of the first matrix.
    /// @tparam DimB Dimension of the second matrix.
    ///
    /// @param A First tridigonal matrix.
    /// @param B Second tridigonal matrix.
    ///
    /// @return Dense tensor-product matrix.
    template<natural_t DimA, natural_t DimB>
    constexpr Matrix<DimA * DimB>
        tensorProduct(
            const tridiagonal_matrix_t<DimA>& A,
            const tridiagonal_matrix_t<DimB>& B) noexcept
    {
        constexpr natural_t ResultDim = DimA * DimB;

        Matrix<ResultDim> Result{};

        // Helper lambda for reading tridigonal elements
        auto getTridigonalElement =
            [](const auto& M,
                natural_t row,
                natural_t col) constexpr -> complex_t
        {
            if (row == col)
            {
                return M[MAINDIAGONAL][row];
            }

            if (col == row + 1)
            {
                return M[SUPERDIAGONAL][row];
            }

            if (row == col + 1)
            {
                return M[SUBDIAGONAL][row];
            }

            return complex_t::zero();
        };

        for (natural_t i = 0; i < DimA; ++i)
        {
            for (natural_t j = 0; j < DimA; ++j)
            {
                const complex_t ContrnA =
                    getTridigonalElement(A, i, j);

                // Skip zero elements
                if (ContrnA == complex_t::zero())
                {
                    continue;
                }

                for (natural_t n = 0; n < DimB; ++n)
                {
                    for (natural_t m = 0; m < DimB; ++m)
                    {
                        const complex_t ContribB =
                            getTridigonalElement(B, n, m);

                        // Skip zero elements
                        if (ContribB == complex_t::zero())
                        {
                            continue;
                        }

                        const natural_t Row = i * DimB + n;
                        const natural_t Col = j * DimB + m;
                        Result.at(Row, Col) = ContrnA * ContribB;
                    }
                }
            }
        }

        return Result;
    }

    template<natural_t Dim>
    constexpr void make_kronecker_sum_indices(
        std::array<std::array<natural_t, 5>, Dim*Dim>& indices,
        std::array<natural_t, Dim*Dim>& counts) noexcept
    {
        constexpr natural_t FullDim = Dim * Dim;

        for (natural_t row = 0; row < FullDim; ++row)
        {
            const natural_t i = row / Dim;
            const natural_t j = row % Dim;

            natural_t count = 0;

            // center
            indices[row][count++] = row;

            // (i-1, j)
            if (i > 0)
            {
                indices[row][count++] = (i - 1) * Dim + j;
            }

            // (i+1, j)
            if (i + 1 < Dim)
            {
                indices[row][count++] = (i + 1) * Dim + j;
            }

            // (i, j-1)
            if (j > 0)
            {
                indices[row][count++] = i * Dim + (j - 1);
            }

            // (i, j+1)
            if (j + 1 < Dim)
            {
                indices[row][count++] = i * Dim + (j + 1);
            }

            counts[row] = count;
        }
    }
}