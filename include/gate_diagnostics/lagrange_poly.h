#pragma once
#include "core_types.h"


namespace KetCat
{

    /// @brief Exact polynomial interpolator using Newton divided differences.
    ///
    /// @details
    ///    Constructs an interpolation polynomial:
    ///
    ///      P(xᵢ) = yᵢ
    ///
    ///    for every supplied sample point.
    ///
    ///    Newton form:
    ///
    ///      P(x) = a₀
    ///           + a₁(x - x₀)
    ///           + a₂(x - x₀)(x - x₁)
    ///           + ...
    ///
    ///    The divided-difference formulation is numerically more stable
    ///    than directly expanding the polynomial basis.
    ///
    ///    Complexity:
    ///      fit()      → O(n²)
    ///      evaluate() → O(n)
    ///
    /// @note
    ///    Interpolation is exact at the provided anchor points.
    ///
    /// @warning
    ///    Very high polynomial orders may exhibit Runge oscillation.
    template <natural_t DataPoints>
    class PolynomialInterpolator
    {
        ///@brief Interpolation anchor x coordinates.
        std::array<real_t, DataPoints> m_X;

        ///@brief Newton divided-difference coefficients.
        std::array<real_t, DataPoints> m_Coefficients;

    public:

        /// @brief Fit interpolation polynomial through sample points.
        ///
        /// @param x Sample x coordinates.
        /// @param y Sample y coordinates.
        constexpr void PolynomialInterpolator(const std::array<real_t, DataPoints>& x, const std::array<real_t, DataPoints>& y) noexcept
            : m_X(x), m_Coefficients(y)
        {
            // Build upper-triangular divided-difference table in-place.
            //
            // After completion:
            //   m_coefficients[i] contains Newton coefficient aᵢ.
            for (natural_t Order = 1; Order < n; ++Order)
            {
                for (natural_t i = DataPoints; i >= Order; --i)
                {
                    const real_t dx =  m_x[i] - m_x[i - Order];

                    m_Coefficients[i] = (m_Coefficients[i] - m_Coefficients[i - 1]) / dx;

                    // Prevent unsigned wraparound.
                    if (i == Order)
                    {
                        break;
                    }
                }
            }
        }

        /// @brief Evaluate interpolated polynomial at arbitrary x.
        ///
        /// @param x Query coordinate.
        ///
        /// @return Interpolated y value.
        ///
        /// @throws std::runtime_error
        ///    If fit() has not been called.
        ///
        /// @details
        ///    Evaluates Newton polynomial using nested multiplication:
        ///
        ///      P(x) =
        ///        a₀
        ///        + (x-x₀)[
        ///            a₁
        ///            + (x-x₁)[
        ///                a₂ + ...
        ///              ]
        ///          ]
        ///
        ///    Equivalent to Horner evaluation in Newton basis.
        real_t evaluate(real_t x) const
        {
             // Start from highest-order coefficient.
            real_t Result = m_Coefficients[DataPoints - 1];

            // Reverse nested multiplication.
            for (natural_t i = DataPoints - 1; i-- > 0;)
            {
                Result = Result * (x - m_X[i]) + m_Coefficients[i];
            }

            return Result;
        }
    };
}
