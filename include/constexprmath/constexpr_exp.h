#pragma once
#include <concepts>
#include "constexpr_complex.h"

namespace ConstexprMath
{
    /// @brief Compute the exponential function using Taylor series expansion.
    /// Now generalized to my own Constexpt Complex type as well.
    /// Moved into a separate file to avoid circular includes (Complex depends on
    /// trigon, trigon depends on core_functions).
    /// @param x The exponent value.
    /// @param N The number of terms in the Taylor series (default is 20).
    template <unsigned int Terms, typename FloatOrComplex>
        requires std::floating_point<FloatOrComplex> ||
    std::same_as<FloatOrComplex, ConstexprMath::Complex<typename FloatOrComplex::FloatType>>
        constexpr FloatOrComplex exp(FloatOrComplex x)
    {
        FloatOrComplex Sum = FloatOrComplex(1.0);
        FloatOrComplex Term = FloatOrComplex(1.0);

        for (unsigned n = 1; n <= Terms; ++n)
        {
            Term = Term * (x / static_cast<double>(n));
            Sum = Sum + Term;
        }

        return Sum;
    }
}