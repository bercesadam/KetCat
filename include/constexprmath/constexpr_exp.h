#pragma once
#include <concepts>
#include "constexpr_complex.h"

namespace ConstexprMath
{
    /// @brief Compute the real exponential function exp(x) for real-valued arguments.
    ///
    /// Uses range reduction followed by a Taylor series expansion to ensure
    /// accurate results for arbitrarily large or small arguments.
    ///
    /// Algorithm:
    ///   1. Split x into integer part n and fractional part frac ∈ [0, 1):
    ///          x = n + frac
    ///   2. Compute exp(frac) via Taylor series (converges fast for frac ∈ [0,1)):
    ///          exp(frac) = Σ_{k=0}^{Terms} frac^k / k!
    ///   3. Compute e^n by repeated multiplication of the Euler number E.
    ///   4. Combine:
    ///          exp(x) = e^n * exp(frac)         if n >= 0
    ///          exp(x) = exp(frac) / e^|n|       if n < 0
    ///
    /// This avoids the poor convergence of a naive Taylor series for |x| >> 1,
    /// which would require an impractically large number of terms.
    ///
    /// @tparam Terms  Number of Taylor series terms for the fractional part.
    ///                20 terms gives error < 1e-17 for frac ∈ [0, 1).
    /// @param x       Real-valued exponent (any finite value).
    /// @return        exp(x) to high precision.
    template <std::floating_point FloatType, unsigned int Terms = 30>
    static constexpr FloatType exp(FloatType x) noexcept
    {
        constexpr FloatType E = 2.718281828459045235360287;

        int n = static_cast<int>(x);
        FloatType frac  = x - static_cast<FloatType>(n);

        // Taylor series for exp(frac), frac ∈ [0, 1)
        FloatType Sum  = 1.0;
        FloatType Term = 1.0;
        for (unsigned k = 1; k <= Terms; ++k)
        {
            Term *= frac / static_cast<FloatType>(k);
            Sum  += Term;
        }

        // e^|n| by repeated multiplication
        int absN = (n >= 0) ? n : -n;
        FloatType ePowN = 1.0;
        for (int k = 0; k < absN; ++k)
            ePowN *= E;

        return (n >= 0) ? ePowN * Sum : Sum / ePowN;
    }

    /// @brief Compute the complex exponential function exp(z) for complex-valued arguments.
    ///
    /// Uses Euler's formula to decompose the complex exponential exactly:
    ///
    ///     exp(a + ib) = exp(a) · (cos(b) + i·sin(b))
    ///
    /// This is both mathematically exact and numerically superior to applying
    /// a Taylor series directly to a complex argument, because:
    ///   - exp(a) is computed via the real overload with range reduction.
    ///   - cos(b) and sin(b) are computed via their own constexpr implementations.
    ///   - No convergence issues arise regardless of the magnitude of b,
    ///     which is important for phase factors e^{-iEt} where E·t can be large.
    ///
    /// @tparam Terms  Number of Taylor series terms forwarded to the real exp overload.
    ///                Only affects precision of the real part magnitude.
    /// @param x       Complex-valued exponent z = a + ib.
    /// @return        exp(z) as a complex number.
    template <std::floating_point FloatType, unsigned int Terms = 30>
    static constexpr ConstexprMath::Complex<FloatType> exp(ConstexprMath::Complex<FloatType> x) noexcept
    {
        FloatType Magnitude = exp<FloatType, Terms>(x.re);
        return ConstexprMath::Complex<FloatType>
        (
            Magnitude * cos(x.im),
            Magnitude * sin(x.im)
        );
    }

    /// @brief Compute the natural logarithm log(x) for x > 0 at compile time.
    ///
    /// Uses range reduction to map x into m ∈ [1, 2), then evaluates:
    ///
    ///     log(x) = k · log(2) + log(m)
    ///
    /// where k is the integer exponent such that x = 2^k · m.
    ///
    /// log(m) is then computed via the atanh identity:
    ///
    ///     log(m) = 2 · atanh((m - 1) / (m + 1))
    ///
    /// using the series:
    ///
    ///     atanh(y) = Σ_{n=0}^{Terms} y^{2n+1} / (2n+1)     |y| < 1
    ///
    /// For m ∈ [1, 2), the argument y = (m-1)/(m+1) ∈ [0, 1/3),
    /// so the series converges geometrically fast — 30 terms gives
    /// error below 1e-15.
    ///
    /// @tparam Terms  Number of atanh series terms.
    ///                30 terms gives error < 1e-15 for all m ∈ [1, 2).
    /// @param x       Real-valued argument, must be > 0.
    /// @return        Natural logarithm of x.
    template <std::floating_point FloatType, unsigned int Terms = 30>
    static constexpr FloatType log(FloatType x) noexcept
    {
        constexpr FloatType Log2 = 0.6931471805599453094172321;

        // Range reduction: find k such that x = 2^k * m, m ∈ [1, 2)
        int k    = 0;
        FloatType m = x;

        while (m >= 2.0) { m *= 0.5; ++k; }
        while (m <  1.0) { m *= 2.0; --k; }

        // atanh series: y = (m-1)/(m+1) ∈ [0, 1/3) for m ∈ [1, 2)
        FloatType y    = (m - 1.0) / (m + 1.0);
        FloatType y2   = y * y;
        FloatType Term = y;
        FloatType Sum  = y;

        for (unsigned n = 1; n <= Terms; ++n)
        {
            Term *= y2;
            Sum  += Term / static_cast<FloatType>(2 * n + 1);
        }

        return static_cast<FloatType>(k) * Log2 + 2.0 * Sum;
    }

    /// @brief Compute base^exponent for real-valued base and exponent.
    ///
    /// Uses the identity:
    ///
    ///     base^exponent = exp(exponent · log(base))
    ///
    /// Both exp() and log() use range reduction internally, so this is
    /// accurate for the full range of values encountered on logarithmic
    /// spatial grids (e.g. base = Extent/RMin ≈ 130000, exponent ∈ [0, 1]).
    ///
    /// @tparam Terms  Number of series terms forwarded to exp() and log().
    /// @param base      Must be > 0.
    /// @param exponent  Any real value.
    /// @return          base^exponent.
    template <std::floating_point FloatType, unsigned int Terms = 30>
    static constexpr FloatType pow(FloatType base, FloatType exponent) noexcept
    {
        return exp<FloatType, Terms>(exponent * log<FloatType, Terms>(base));
    }
}