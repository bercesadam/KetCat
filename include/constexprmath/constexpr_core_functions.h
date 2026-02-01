#pragma once
#include <concepts>
#include <limits>
#include <cmath>


/// @file
/// @brief Small constexpr integer utilities used for dimensions and bit manipulations.
///
/**
 * @details
 * This header provides:
 *  - `pow2(n)` : compute 2^n at compile time via bit-shift (returns an unsigned integral).
 *  - `is_power_of_two(x)` : test whether an integer is a (positive) power of two.
 *
 * Both functions are constexpr and noexcept and intended for use in compile-time
 * dimension computations (matrix sizes, state vector lengths) and static assertions.
 */

namespace ConstexprMath
{
    /**
     * @brief Compute 2 raised to the power of `n` at compile time.
     *
     * @tparam UIntType  Unsigned integral type used for the exponent and result.
     * @param n          Exponent (non-negative).
     * @return 2^n as `UIntType`.
     */
    template <std::unsigned_integral UIntType>
    constexpr UIntType pow2(UIntType n) noexcept
    {
        // Left-shift 1 by n positions: 1 << n == 2^n for unsigned types.
        return UIntType{ 1 } << n;
    }

    /**
     * @brief Determine whether a value is an exact power of two.
     *
     * @tparam UIntType  Unsigned integral type for the argument.
     * @param x          Value to test.
     * @return true if x is a power of two (1, 2, 4, 8, ...); false otherwise.
     *
     * @note This uses the classic bit trick: x > 0 && (x & (x - 1)) == 0.
     *       The operation is constexpr and does not allocate memory.
     */
    template <std::unsigned_integral UIntType>
    constexpr bool isPowerOfTwo(UIntType x) noexcept
    {
        // Fast test: powers of two have exactly one bit set.
        return x > 0 && (x & (x - 1)) == 0;
    }

    /**
     * @brief Compute the square root of a floating-point number at compile time.
     *
     * @param x  The input value (must be non-negative).
     * @return   The square root of x, or NaN if x is negative.
     *
     * @note This implementation uses the Newton-Raphson method for computing
     *       the square root and is constexpr-friendly.
	 */
    template <std::floating_point FloatType>
    constexpr FloatType sqrt(FloatType x)
    {
        // Handle edge cases
        if (x < 0.0) return std::numeric_limits<FloatType>::quiet_NaN();
        if (x == 0.0 || x == std::numeric_limits<FloatType>::infinity()) return x;

        // Recursive lambda for Newton-Raphson
        auto sqrtRec = [](FloatType x, FloatType curr, FloatType prev, auto&& self) -> FloatType
        {
            return curr == prev ? curr
                : self(x, 0.5 * (curr + x / curr), curr, self);
        };

        return sqrtRec(x, x, 0.0, sqrtRec);
    }

    /// @brief Compute the exponential function using Taylor series expansion.
    /// @param x The exponent value.
    /// @param N The number of terms in the Taylor series (default is 20).
    template <unsigned int Terms, std::floating_point FloatType>
    constexpr FloatType exp(FloatType x)
    {
        FloatType sum = 1.0;
        FloatType term = 1.0;
        for (unsigned n = 1; n <= Terms; ++n)
        {
            term *= x / n;
            sum += term;
        }
        return sum;
    }

    /// @brief Compute absolute value of a number
    /// @param x The input value 
    template <std::floating_point FloatType>
    constexpr FloatType abs(FloatType x)
    {
        return (x < 0.0) ? -x : x;
    };

	/// @brief Compute the factorial of a non-negative integer at compile time.
	/// @tparam IntType  Integral type for the argument and result.
	/// @param n         Non-negative integer.
	template <std::unsigned_integral UIntType>
    constexpr UIntType factorial(UIntType n) noexcept
    {
        return (n <= 1) ? 1 : n * factorial(n - 1);
	}
}