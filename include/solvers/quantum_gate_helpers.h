#pragma once
#include "core_types.h"

namespace KetCat::QCC
{
    /// @file
    /// @brief Helpers for applying gate matrices and validating matrix properties.
    /// @details
    /// Provides:
    ///  - Type trait `is_gate_matrix` to identify gate matrix types.
    ///  - `apply_unitary`: constexpr matrix–vector multiplication.
    ///  - `is_valid_square_matrix`: compile‑time check for square matrices of size 2^k.
    ///  - `is_unitary`: constexpr‑checkable verification that a matrix is unitary.

    /// @brief Type trait to check if a type is a gate matrix — a square array of complex
    ///        numbers whose size (N × N) is a power of two.
    template<typename T>
    struct is_gate_matrix : std::false_type {};

    /// @brief Specialization for square matrices represented as nested std::array.
    template<dimension_t N>
    struct is_gate_matrix<std::array<std::array<cplx_t, N>, N>>
    {
        static constexpr dimension_t dim = N;
        static constexpr bool value = ConstexprMath::isPowerOfTwo(N);
    };

    /// @brief Convenience `_v` alias for the `is_gate_matrix` trait.
    template<typename T>
    inline constexpr bool is_gate_matrix_v = is_gate_matrix<T>::value;

    /// @brief Check whether a provided square complex matrix is unitary.
    /// @tparam M A constexpr matrix value (template non‑type parameter).
    /// @return True if M · M† ≈ I, false otherwise.
    /// @details
    /// Computes the Hermitian adjoint M† of M, multiplies M by M†, and verifies that the
    /// result matches the identity matrix within a small tolerance.  
    /// This function is intended for compile‑time validation of gate matrices, and
    /// performs exact equality checks on real and imaginary components except for an
    /// epsilon tolerance.
    template<auto M>
    constexpr bool is_unitary()
    {
        // Matrix dimension from the type trait.
        constexpr dimension_t Dim =
            is_gate_matrix<std::remove_cvref_t<decltype(M)>>::dim;

        // Tolerance for floating‑point comparison.
        constexpr real_t Epsilon = 1E-9;

        // Check unitarity: M · M† == I.
        for (dimension_t i = 0; i < Dim; i++)
            for (dimension_t j = 0; j < Dim; j++)
            {
                cplx_t Sum{ 0,0 };

                for (dimension_t k = 0; k < Dim; k++)
                    Sum = Sum + M[k][i].conj() * M[k][j];

                if (i == j &&
                    (ConstexprMath::abs(Sum.re - 1) > Epsilon || ConstexprMath::abs(Sum.im) > Epsilon))
                    return false;

                if (i != j &&
                    (ConstexprMath::abs(Sum.re) > Epsilon || ConstexprMath::abs(Sum.im) > Epsilon))
                    return false;
            }

        return true;
    }
}
