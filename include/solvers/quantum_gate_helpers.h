#pragma once
#include "core_types.h"

namespace KetCat::QCC
{
    /// @file
    /// @brief Helpers for applying gate matrices and validating matrix properties.
    ///
    /**
     * @details
     * This header provides:
     *  - Type trait `is_gate_matrix` to identify gate matrix types.
     *  - `apply_unitary` : constexpr matrix-vector multiplication (used to apply a local gate).
     *  - `is_valid_square_matrix` : compile-time check for square matrices with power-of-two size.
     *  - `is_unitary` : constexpr runtime/checkable check that a matrix is unitary.
     */

     /// @brief  Type trait to check if a type is a gate matrix:
     /// square array of complex numbers where the dimension is a power of two.
    template<typename T>
    struct is_gate_matrix : std::false_type {};

    // Specialization for square matrices represented as nested std::array
    template<dimension_t N>
    struct is_gate_matrix<std::array<std::array<cplx_t, N>, N>> {
        static constexpr dimension_t dim = N;
        static constexpr bool value = ConstexprMath::isPowerOfTwo(N);
    };

    /// C++17 style '_v' helper
    template<typename T>
    inline constexpr bool is_gate_matrix_v = is_gate_matrix<T>::value;

    /// @brief  Checks whether a provided square complex matrix is unitary.
    /// @tparam Dim     Dimension of the square matrix (2^k).
    /// @param mat      The matrix to check.
    /// @return         True if mat * mat^� equals the identity, false otherwise.
    ///
    /// @details
    /// This function computes the conjugate transpose (Hermitian adjoint) of `mat`,
    /// multiplies `mat` by its conjugate transpose and verifies that the product is
    /// equal to the identity matrix within exact arithmetic. Because this routine
    /// uses exact float_t checks for equality of re/im components, it is intended
    /// for constexpr / compile-time constructed matrices used as gates.
    template<auto M>
    constexpr bool is_unitary()
    {
        // Get the matrix dimension from the type trait
        constexpr dimension_t Dim = is_gate_matrix<std::remove_cvref_t<decltype(M)>>::dim;

        // Epsilon for floating-point comparison
        constexpr float_t Epsilon = 1E-9;

        // Check unitarity: M * M^� == I
        for (dimension_t i = 0; i < Dim; i++)
            for (dimension_t j = 0; j < Dim; j++)
            {
                cplx_t Sum{ 0,0 };
                for (dimension_t k = 0; k < Dim; k++)
                    Sum = Sum + M[k][i].conj() * M[k][j];

                constexpr auto abs = [](float_t v) -> float_t
                {
                    return (v < 0.0) ? -v : v;
                };

                if (i == j && (abs(Sum.re - 1) > Epsilon || abs(Sum.im) > Epsilon))
                    return false;
                if (i != j && (abs(Sum.re) > Epsilon || abs(Sum.im) > Epsilon))
                    return false;
            }
        return true;
    }
}