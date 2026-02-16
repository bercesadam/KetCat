#pragma once
#include <cstdint>
#include "constexpr_core_functions.h"

/// @file
/// @brief Constexpr-friendly trigonometric helpers used to build unitary gates.
/// @details
/// This header provides constexpr implementations of sine, cosine, arctangent,
/// and related helpers needed when constructing compile‑time quantum gate
/// matrices (e.g., QFT / IQFT).
///
/// Key features:
/// 1. Proper quadrant-aware range reduction so sin(x) and cos(x) have correct signs.
/// 2. Higher‑order Taylor‑style polynomials for improved accuracy on [-π/4, π/4].
/// 3. Fully constexpr‑friendly — no dynamic allocations, suitable for
///    compile‑time matrix construction.
///
/// Polynomial notes:
///   sin(x) ≈ x - x³/3! + x⁵/5! - x⁷/7! + x⁹/9! - x¹¹/11!
///   cos(x) ≈ 1 - x²/2! + x⁴/4! - x⁶/6! + x⁸/8! - x¹⁰/10! + x¹²/12!
/// Range reduction maps angles into [-π/4, π/4] with quadrant handling.
/// Works for any double input in a constexpr context.

namespace ConstexprMath
{
    /// @brief Mathematical constant π in double precision.
    static constexpr double Pi = 3.141592653589793238462643383279502884;

    /// @brief Constexpr factorial (used for defining polynomial coefficients).
    /// @param n Non-negative integer.
    /// @return n! as a double.
    constexpr double factorial(int n) noexcept
    {
        double r = 1.0;
        for (int i = 2; i <= n; ++i) r *= i;
        return r;
    }

    /// @brief Constexpr floor implementation for double.
    /// @param x Input value.
    /// @return The largest integer ≤ x.
    constexpr int floorConstexpr(double x) noexcept
    {
        int32_t i = static_cast<int32_t>(x);
        return (x < static_cast<double>(i)) ? (i - 1) : i;
    }

    /// @brief Polynomial approximation to sin(x) around 0.
    /// @details Accurate for |x| ≤ π/4.
    constexpr double sinPoly(double x) noexcept
    {
        const double x2 = x * x;
        double result = x;
        double term = x;
        term *= -x2 / 6.0;  result += term;    // -x^3/3!
        term *= -x2 / 20.0; result += term;    // +x^5/5!
        term *= -x2 / 42.0; result += term;    // -x^7/7!
        term *= -x2 / 72.0; result += term;    // +x^9/9!
        term *= -x2 / 110.0; result += term;   // -x^11/11!
        return result;
    }

    /// @brief Polynomial approximation to cos(x) around 0.
    /// @details Accurate for |x| ≤ π/4.
    constexpr double cosPoly(double x) noexcept
    {
        const double x2 = x * x;
        double result = 1.0;
        double term = -x2 / 2.0;  result += term;    // -x^2/2!
        term *= -x2 / 12.0; result += term;          // +x^4/4!
        term *= -x2 / 30.0; result += term;          // -x^6/6!
        term *= -x2 / 56.0; result += term;          // +x^8/8!
        term *= -x2 / 90.0; result += term;          // -x^10/10!
        term *= -x2 / 132.0; result += term;         // +x^12/12!
        return result;
    }

    /// @brief Reduce angle to [-π/4, π/4] and determine quadrant.
    /// @param x Input angle in radians.
    /// @param xr Output reduced angle.
    /// @return Quadrant index in [0, 3].
    constexpr int reduceQuadrant(double x, double& xr) noexcept
    {
        // Map x = k*(π/2) + xr, with xr in [-π/4, π/4].
        double k = floorConstexpr((x + Pi / 4.0) / (Pi / 2.0));
        xr = x - k * (Pi / 2.0);
        return static_cast<int>(k) & 3; // quadrant mod 4
    }

    /// @brief Constexpr sine with quadrant-aware range reduction.
    /// @param x Angle in radians.
    constexpr double sin(double x) noexcept
    {
        double xr;
        int q = reduceQuadrant(x, xr);
        switch (q)
        {
        case 0: return  sinPoly(xr);
        case 1: return  cosPoly(xr);
        case 2: return -sinPoly(xr);
        case 3: return -cosPoly(xr);
        }
        return 0.0; // unreachable
    }

    /// @brief Constexpr cosine with quadrant-aware range reduction.
    /// @param x Angle in radians.
    constexpr double cos(double x) noexcept
    {
        double xr;
        int q = reduceQuadrant(x, xr);
        switch (q)
        {
        case 0: return  cosPoly(xr);
        case 1: return -sinPoly(xr);
        case 2: return -cosPoly(xr);
        case 3: return  sinPoly(xr);
        }
        return 0.0; // unreachable
    }

    /// @brief Polynomial approximation for atan(z) for |z| ≤ 1.
    constexpr double atanPoly(double z) noexcept
    {
        double z2 = z * z;
        double term = z;
        double result = z;

        term *= -z2 / 3.0;        result += term;  // -z^3/3
        term *= -z2 * 3.0 / 5.0;  result += term;  // +z^5/5
        term *= -z2 * 5.0 / 7.0;  result += term;  // -z^7/7
        term *= -z2 * 7.0 / 9.0;  result += term;  // +z^9/9
        return result;
    }

    /// @brief Constexpr arctangent with range reduction.
    constexpr double atan(double z) noexcept
    {
        if (z > 1.0)  return  Pi / 2.0 - atanPoly(1.0 / z);
        if (z < -1.0) return -Pi / 2.0 - atanPoly(1.0 / z);
        return atanPoly(z);
    }

    /// @brief Constexpr atan2 implementation for all quadrants.
    constexpr double atan2(double y, double x) noexcept
    {
        if (x > 0.0) return atan(y / x);
        if (x < 0.0)
        {
            return (y >= 0.0)
                ? atan(y / x) + Pi
                : atan(y / x) - Pi;
        }
        // x == 0
        if (y > 0.0) return  Pi / 2.0;
        if (y < 0.0) return -Pi / 2.0;
        return 0.0; // undefined
    }

    /// @brief Constexpr arcsine using atan2 formulation.
    /// @param x Input value clamped to [-1, 1].
    /// @return asin(x) in radians.
    constexpr double asin(double x) noexcept
    {
        if (x <= -1.0) return -Pi / 2.0;   // asin(-1) = -π/2
        if (x >= 1.0)  return  Pi / 2.0;   // asin( 1) =  π/2

        // asin(x) = atan2(x, sqrt(1 - x^2))
        return atan2(x, sqrt(1.0 - x * x));
    }

    /// @brief Constexpr arccos using atan2 formulation.
    /// @param x Input value clamped to [-1, 1].
    constexpr double acos(double x) noexcept
    {
        if (x <= -1.0) return Pi;
        if (x >= 1.0)  return 0.0;
        return atan2(sqrt(1.0 - x * x), x);
    }
}
