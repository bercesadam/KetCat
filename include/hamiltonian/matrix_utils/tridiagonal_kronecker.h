#pragma once
#include "core_types.h"
#include 

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
    constexpr square_matrix_t<DimA* DimB>
        tensorProduct(
            const tridigonal_matrix_t<DimA>& A,
            const tridigonal_matrix_t<DimB>& B) noexcept
    {
        constexpr natural_t ResultDim = DimA * DimB;

        square_matrix_t<ResultDim> Result{};

        // Helper lambda for reading tridigonal elements
        auto getTridigonalElement =
            [](const auto& M,
                natural_t row,
                natural_t col) constexpr -> complex_t
        {
            if (row == col)
            {
                return M[MAINDiGONAL][row];
            }

            if (col == row + 1)
            {
                return M[SUPERDiGONAL][row];
            }

            if (row == col + 1)
            {
                return M[SUBDiGONAL][row];
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
                        Result[Row][Col] = ContrnA * ContribB;
                    }
                }
            }
        }

        return Result;
    }
}