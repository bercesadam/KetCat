#pragma once
#include <concepts>
#include <limits>
#include <cmath>

/// @file
/// @brief Small constexpr integer utilities used for dimensions and bit manipulations.
/// @details
/// Provides several small constexpr utilities:
///   - `pow2(n)`: compute 2^n at compile time using bit shift.
///   - `isPowerOfTwo(x)`: test whether an integer is an exact power of two.
///   - `sqrt(x)`: constexpr Newton–Raphson square root.
///   - `abs(x)`: constexpr absolute value.
///   - `factorial(n)`: constexpr factorial for unsigned integrals.
/// Intended for dimension computations, static assertions, and other compile‑time operations.

namespace ConstexprMath
{
    /// @brief Compute 2 raised to the power of `n` at compile time.
    /// @tparam UIntType Unsigned integral type.
    /// @param n Exponent (non‑negative).
    /// @return 2^n as `UIntType`.
    template <std::unsigned_integral UIntType>
    constexpr UIntType pow2(UIntType n) noexcept
    {
        // Left-shift 1 by n bits: 1 << n == 2^n for unsigned types.
        return UIntType{ 1 } << n;
    }

    
    // @brief Compute `base` raised to the power of `exp` at compile time using
    ///        exponentiation by squaring.
    /// @tparam T Arithmetic type (integer or floating‑point).
    /// @param base The base value.
    /// @param exp  Exponent (non‑negative).
    /// @return base^exp as type `T`.
    template <std::unsigned_integral UIntType>
    constexpr UIntType pow(UIntType n, UIntType exponent) noexcept
    {
        UIntType Result = 1;
        while (exponent > 0)
        {
            if (exponent & 1)
            {
                Result *= n;
            }
            n *= n;
            exponent >>= 1;
        }
        return Result;
    }

    /// @brief Determine whether a value is an exact power of two.
    /// @tparam UIntType Unsigned integral type.
    /// @param x Value to check.
    /// @return true if x is a power of two (1,2,4,...); false otherwise.
    /// @note Uses the classic bit trick: x > 0 && (x & (x - 1)) == 0.
    template <std::unsigned_integral UIntType>
    constexpr bool isPowerOfTwo(UIntType x) noexcept
    {
        // Powers of two have exactly one bit set.
        return x > 0 && (x & (x - 1)) == 0;
    }

    /// @brief Compute square root at compile time using Newton–Raphson iteration.
    /// @tparam FloatType Floating‑point type.
    /// @param x Input value (must be non‑negative).
    /// @return sqrt(x), or NaN if x is negative.
    template <std::floating_point FloatType>
    constexpr FloatType sqrt(FloatType x)
    {
        // Handle edge cases.
        if (x < 0.0) return std::numeric_limits<FloatType>::quiet_NaN();
        if (x == 0.0 || x == std::numeric_limits<FloatType>::infinity()) return x;

        // Recursive Newton–Raphson lambda.
        auto sqrtRec = [](FloatType x, FloatType curr, FloatType prev, auto&& self) -> FloatType
        {
            return curr == prev
                ? curr
                : self(x, 0.5 * (curr + x / curr), curr, self);
        };

        return sqrtRec(x, x, 0.0, sqrtRec);
    }

    /// @brief Compute absolute value of a number.
    /// @tparam NumericType Arithmetic type.
    /// @param x Input value.
    template <typename NumericType>
        requires std::is_arithmetic_v<NumericType>
    constexpr NumericType abs(NumericType x)
    {
        return (x < NumericType{}) ? -x : x;
    }

    /// @brief Compute factorial of a non‑negative integer at compile time.
    /// @tparam UIntType Unsigned integral type.
    /// @param n Non‑negative integer.
    /// @return n! computed recursively at compile time.
    template <std::unsigned_integral UIntType>
    constexpr UIntType factorial(UIntType n) noexcept
    {
        return (n <= 1) ? 1 : n * factorial(n - 1);
    }
}
