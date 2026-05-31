#pragma once
#include "core_types.h"
#include "backend_traits.h"
#include <type_traits>


namespace KetCat
{
    // =========================================================================
    // Matrix LevelCountension Extractor Traits
    // =========================================================================

    /// @brief Primary template to extract dimension from banded matrix types.
    template<typename T>
    struct MatrixTraits;

    /// @brief Specialization to extract LevelCount from tridiagonal_matrix_t.
    template<natural_t LevelCount>
    struct MatrixTraits<tridiagonal_matrix_t<LevelCount>>
    {
        static constexpr natural_t LevelCountVal = LevelCount;
    };

    /// @brief Specialization to extract LevelCount from five_band_matrix_t.
    template<natural_t LevelCount>
    struct MatrixTraits<five_band_matrix_t<LevelCount>>
    {
        static constexpr natural_t LevelCountVal = LevelCount * LevelCount;
    };

    // =========================================================================
    // Matrix Access Interface
    // =========================================================================

    /// @brief Matrix access traits to provide a uniform interface for accessing elements of different matrix types.
    template<typename MatrixType>
    struct MatrixAccess;

    // =========================================================================
    // Specialization for TRIDIAGONAL Matrix
    // =========================================================================
    template<natural_t LevelCount>
    struct MatrixAccess<tridiagonal_matrix_t<LevelCount>>
    {
        /// @brief Unified element getter for tridiagonal structures.
        static constexpr complex_t get(
            const tridiagonal_matrix_t<LevelCount>& M,
            natural_t row, natural_t col) noexcept
        {
            if (row == col)     return M[MAINDIAGONAL][row];
            if (col == row + 1) return M[SUPERDIAGONAL][row];
            if (row == col + 1) return M[SUBDIAGONAL][row];

            return complex_t::zero();
        }

        /// @brief Unified element setter for tridiagonal structures.
        static constexpr void set(
            tridiagonal_matrix_t<LevelCount>& M,
            natural_t row, natural_t col,
            const complex_t& value) noexcept
        {
            if (row == col)     M[MAINDIAGONAL][row] = value;
            if (col == row + 1) M[SUPERDIAGONAL][row] = value;
            if (row == col + 1) M[SUBDIAGONAL][row] = value;
        }
    };

    // =========================================================================
    // Specialization for PENTADIAGONAL Matrix
    // =========================================================================
    template<natural_t LevelCount>
    struct MatrixAccess<five_band_matrix_t<LevelCount>>
    {
        static constexpr complex_t get(const five_band_matrix_t<LevelCount>& M, natural_t row, natural_t col) noexcept {
            if (row == col) return M[MAINDIAGONAL][row];

            // n (belső koordináta) = row % LevelCount
            natural_t n = row % LevelCount;

            if (col == row + 1 && n < LevelCount - 1) return M[SUPERDIAGONAL][row];
            if (row == col + 1 && n > 0)               return M[SUBDIAGONAL][row];
            if (col == row + LevelCount)               return M[UPPER_FAR][row];
            if (row == col + LevelCount)               return M[LOWER_FAR][row];

            return complex_t::zero();
        }

        static constexpr void set(five_band_matrix_t<LevelCount>& M, natural_t row, natural_t col, const complex_t& value) noexcept {
            if (row == col) { M[MAINDIAGONAL][row] = value; return; }

            natural_t n = row % LevelCount;

            if (col == row + 1 && n < LevelCount - 1) { M[SUPERDIAGONAL][row] = value; return; }
            if (row == col + 1 && n > 0) { M[SUBDIAGONAL][row] = value; return; }
            if (col == row + LevelCount) { M[UPPER_FAR][row] = value; return; }
            if (row == col + LevelCount) { M[LOWER_FAR][row] = value; return; }
        }
    };
}
