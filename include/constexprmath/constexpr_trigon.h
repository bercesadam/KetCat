#pragma once
#include <cstdint>
#include "constexpr_core_functions.h"


/// @file
/// @brief Constexpr-friendly trigonometric helpers used to build unitary gates.
///
/// @file
/// @brief Constexpr-friendly sine and cosine with quadrant-aware range reduction
///
/// @details
/// This header provides `constexpr` sine and cosine functions that can be used
/// at compile-time to construct unitary gate matrices like QFT/IQFT.
/// It improves over naive low-order polynomial approximations by:
/// 1) Applying proper quadrant handling so sin(x) and cos(x) have correct signs.
/// 2) Using higher-order Taylor-like polynomials for better accuracy over [-π/4, π/4].
/// 3) Remaining fully constexpr-friendly without any runtime allocations.
///
/// Notes:
///  - The polynomials are truncated Taylor series:
///      sin(x) ≈ x - x^3/6 + x^5/120 - x^7/5040 + x^9/362880 - x^11/39916800
///      cos(x) ≈ 1 - x^2/2 + x^4/24 - x^6/720 + x^8/40320 - x^10/3628800 + x^12/479001600
///  - Range reduction maps any angle to [-π/4, π/4] and handles quadrant permutations.
///  - Works for any double input in a constexpr context.

namespace ConstexprMath
{
    /// @brief Mathematical constant π in double precision
    static constexpr double Pi = 3.141592653589793238462643383279502884;

    /// @brief Factorial function constexpr (optional, used to define coefficients)
    constexpr double factorial(int n) noexcept
    {
        double r = 1.0;
        for (int i = 2; i <= n; ++i) r *= i;
        return r;
    }

    /// @brief Floor function usable in constexpr
    constexpr int floorConstexpr(double x) noexcept
    {
        int32_t i = static_cast<int32_t>(x);
        return (x < static_cast<double>(i)) ? (i - 1) : i;
    }

    /// @brief Polynomial approximation to sin(x) around 0
    /// Accurate for |x| <= π/4
    constexpr double sinPoly(double x) noexcept
    {
        const double x2 = x * x;
        double result = x;
        double term = x;
        term *= -x2 / 6.0; result += term;         // -x^3/3!
        term *= -x2 / 20.0; result += term;        // +x^5/5!
        term *= -x2 / 42.0; result += term;        // -x^7/7!
        term *= -x2 / 72.0; result += term;        // +x^9/9!
        term *= -x2 / 110.0; result += term;       // -x^11/11!
        return result;
    }

    /// @brief Polynomial approximation to cos(x) around 0
    /// Accurate for |x| <= π/4
    constexpr double cosPoly(double x) noexcept
    {
        const double x2 = x * x;
        double result = 1.0;
        double term = -x2 / 2.0; result += term;      // -x^2/2!
        term *= -x2 / 12.0; result += term;           // +x^4/4!
        term *= -x2 / 30.0; result += term;           // -x^6/6!
        term *= -x2 / 56.0; result += term;           // +x^8/8!
        term *= -x2 / 90.0; result += term;           // -x^10/10!
        term *= -x2 / 132.0; result += term;          // +x^12/12!
        return result;
    }

    /// @brief Reduce angle to [-π/4, π/4] and return quadrant index
    /// @param x Input angle in radians
    /// @param xr Output reduced angle (by reference)
    /// @return Quadrant index 0..3
    constexpr int reduceQuadrant(double x, double& xr) noexcept
    {
        // Map x to k*(π/2) + xr, xr in [-π/4, π/4]
        double k = floorConstexpr((x + Pi / 4.0) / (Pi / 2.0));
        xr = x - k * (Pi / 2.0);
        int kmod4 = static_cast<int>(k) & 3; // quadrant modulo 4
        return kmod4;
    }

    /// @brief Constexpr sine with quadrant-aware range reduction
    constexpr double sin(double x) noexcept
    {
        double xr;
        int q = reduceQuadrant(x, xr);
        switch (q)
        {
        case 0: return sinPoly(xr);
        case 1: return cosPoly(xr);
        case 2: return -sinPoly(xr);
        case 3: return -cosPoly(xr);
        }
        return 0.0; // unreachable
    }

    /// @brief Constexpr cosine with quadrant-aware range reduction
    constexpr double cos(double x) noexcept
    {
        double xr;
        int q = reduceQuadrant(x, xr);
        switch (q)
        {
        case 0: return cosPoly(xr);
        case 1: return -sinPoly(xr);
        case 2: return -cosPoly(xr);
        case 3: return sinPoly(xr);
        }
        return 0.0; // unreachable
    }

    // Range-reduced arctangent series (|z| <= 1)
    constexpr double atanPoly(double z) noexcept
    {
        double z2 = z * z;
        double term = z;
        double result = z;

        term *= -z2 / 3.0; result += term;        // -z^3/3
        term *= -z2 * 3.0 / 5.0; result += term;  // +z^5/5
        term *= -z2 * 5.0 / 7.0; result += term;  // -z^7/7
        term *= -z2 * 7.0 / 9.0; result += term;  // +z^9/9
        return result;
    }

    // Full atan with range reduction for |z| > 1
    constexpr double atan(double z) noexcept
    {
        if (z > 1.0)
            return Pi / 2.0 - atanPoly(1.0 / z);
        if (z < -1.0)
            return -Pi / 2.0 - atanPoly(1.0 / z);
        return atanPoly(z);
    }

    // Constexpr atan2 handling all quadrants
    constexpr double atan2(double y, double x) noexcept
    {
        if (x > 0.0) return atan(y / x);
        if (x < 0.0)
        {
            if (y >= 0.0) return atan(y / x) + Pi;
            else return atan(y / x) - Pi;
        }
        // x == 0
        if (y > 0.0) return Pi / 2.0;
        if (y < 0.0) return -Pi / 2.0;
        return 0.0; // undefined, return 0
    }

    constexpr double acos(double x) noexcept
    {
        // Clamp x to [-1,1] to avoid sqrt of negative due to rounding
        if (x <= -1.0) return Pi;
        if (x >= 1.0) return 0.0;
        return atan2(sqrt(1.0 - x * x), x);
    }
}