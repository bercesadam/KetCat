#pragma once
#include "core_types.h"
#include "backend_traits.h"
#include <type_traits>


namespace KetCat
{
    // =========================================================================
    // Matrix Dimension Extractor Traits
    // =========================================================================

    /// @brief Primary template to extract dimension from banded matrix types.
    template<typename T>
    struct MatrixTraits;

    /// @brief Specialization to extract Dim from tridiagonal_matrix_t.
    template<natural_t Dim>
    struct MatrixTraits<tridiagonal_matrix_t<Dim>>
    {
        static constexpr natural_t DimVal = Dim;
    };

    /// @brief Specialization to extract Dim from pentadiagonal_matrix_t.
    template<natural_t Dim>
    struct MatrixTraits<pentadiagonal_matrix_t<Dim>>
    {
        static constexpr natural_t DimVal = Dim;
    };

    // =========================================================================
    // Matrix Concept Definition
    // =========================================================================

    /// @brief Concept to constrain matrix types to supported banded storage representations.
    template<typename T>
    concept any_matrix_t = requires {
        { MatrixTraits<T>::DimVal } -> std::same_as<const natural_t&>;
    } && (
        std::same_as<T, tridiagonal_matrix_t<MatrixTraits<T>::DimVal>> ||
        std::same_as<T, pentadiagonal_matrix_t<MatrixTraits<T>::DimVal>>
        );

    // =========================================================================
    // Matrix Access Interface
    // =========================================================================

    /// @brief Matrix access traits to provide a uniform interface for accessing elements of different matrix types.
    template<typename MatrixType>
    struct MatrixAccess;

    // =========================================================================
    // Specialization for TRIDIAGONAL Matrix
    // =========================================================================
    template<natural_t Dim>
    struct MatrixAccess<tridiagonal_matrix_t<Dim>>
    {
        /// @brief Unified element getter for tridiagonal structures.
        static constexpr complex_t get(
            const tridiagonal_matrix_t<Dim>& M,
            natural_t row, natural_t col) noexcept
        {
            if (row == col)     return M[MAINDIAGONAL][row];
            if (col == row + 1) return M[SUPERDIAGONAL][row];
            if (row == col + 1) return M[SUBDIAGONAL][row];

            return complex_t::zero();
        }

        /// @brief Unified element setter for tridiagonal structures.
        static constexpr void set(
            tridiagonal_matrix_t<Dim>& M,
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
    template<natural_t Dim>
    struct MatrixAccess<pentadiagonal_matrix_t<Dim>>
    {
        /// @brief Unified element getter for pentadiagonal structures.
        static constexpr complex_t get(
            const pentadiagonal_matrix_t<Dim>& M,
            natural_t row, natural_t col) noexcept
        {
            if (row == col)     return M[MAINDIAGONAL][row];
            if (col == row + 1) return M[SUPERDIAGONAL][row];
            if (col == row + 2) return M[SUPERDIAGONAL2][row];
            if (row == col + 1) return M[SUBDIAGONAL][row];
            if (row == col + 2) return M[SUBDIAGONAL2][row];

            return complex_t::zero();
        }

        /// @brief Unified element setter for pentadiagonal structures.
        static constexpr void set(
            pentadiagonal_matrix_t<Dim>& M,
            natural_t row, natural_t col,
            const complex_t& value) noexcept
        {
            if (row == col)     M[MAINDIAGONAL][row] = value;
            if (col == row + 1) M[SUPERDIAGONAL][row] = value;
            if (col == row + 2) M[SUPERDIAGONAL2][row] = value;
            if (row == col + 1) M[SUBDIAGONAL][row] = value;
            if (row == col + 2) M[SUBDIAGONAL2][row] = value;
        }
    };
}
