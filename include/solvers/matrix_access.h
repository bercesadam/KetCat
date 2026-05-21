#pragma once
#include "core_types.h"


namespace KetCat
{
	/// @brief Concept to constrain matrix types to either square or tridiagonal matrices of the same dimension.
    template<typename T>
    concept any_matrix_t =
        std::same_as<T, tridiagonal_matrix_t<T::size()>> ||
        std::same_as<T, square_matrix_t<T::size()>>;

	/// @brief Matrix access traits to provide a uniform interface for accessing elements of different matrix types.
    template<typename MatrixType>
    struct MatrixAccess;

    /// @brief Specialization for square matrices, providing direct access to elements.
    template<natural_t Dim>
    struct MatrixAccess<square_matrix_t<Dim>>
    {
        static constexpr complex_t get (
            const square_matrix_t<Dim>& M,
            natural_t row,
            natural_t col)
        {
            return M[row][col];
        }

        static constexpr void set(
            square_matrix_t<Dim>& M,
            natural_t row,
            natural_t col,
            complex_t value)
        {
            M[row][col] = value;
        }
    };

    template<natural_t Dim>
    struct MatrixAccess<tridiagonal_matrix_t<Dim>>
    {
        static constexpr complex_t get(
            const tridiagonal_matrix_t<Dim>& M,
            natural_t row,
            natural_t col)
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
        }

        static constexpr void set(
            tridiagonal_matrix_t<Dim>& M,
            natural_t row,
            natural_t col,
            complex_t value)
        {
            if (row == col)
            {
                M[MAINDIAGONAL][row] = value;
            }

            else if (col == row + 1)
            {
                M[SUPERDIAGONAL][row] = value;
            }

            else if (row == col + 1)
            {
                M[SUBDIAGONAL][row] = value;
            }
        }
    };
}
