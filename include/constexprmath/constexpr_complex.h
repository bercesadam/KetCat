#pragma once
#include <concepts>
#include "constexpr_trigon.h"


/// @file
/// @brief Constexpr-capable complex number type used throughout the project.
///
/**
 * @details
 * This header defines `ConstexprMath::Complex<T>`, a small constexpr-friendly
 * complex number type with a minimal set of operations required by the
 * simulator (construction, addition, subtraction, multiplication, in-place
 * accumulation and conjugation). All members and functions are `constexpr`
 * and `noexcept` where appropriate so they can be used in compile-time
 * contexts (constexpr initialization of gate matrices, etc.).
 */

namespace ConstexprMath
{
    /**
     * @brief Simple constexpr complex number for numeric computations.
     *
     * @tparam T  Floating-point type to store real and imaginary parts (e.g. real_t).
     *
     * The type exposes public data members `re` and `im` and lightweight
     * arithmetic operators implemented in a constexpr-friendly way. It is
     * intentionally minimal (no exceptions, no heap allocations).
     */
    template <std::floating_point _FloatType>
    struct Complex
    {
        using FloatType = _FloatType;

        /// @brief The real part of the complex number.
        FloatType re = 0.0;
        /// @brief The imaginary part of the complex number.
        FloatType im = 0.0;

        /// @brief Default constructor producing 0 + 0i.
        constexpr Complex() noexcept = default;

        /// @brief Construct a complex number from a real value (imaginary part = 0).
        /// @param r  Real part.
        constexpr Complex(FloatType r) noexcept : re(r), im(0.0) {}

        /// @brief Construct a complex number from real and imaginary parts.
        /// @param r  Real part.
        /// @param i  Imaginary part.
        constexpr Complex(FloatType r, FloatType i) noexcept : re(r), im(i) {}

        /// @brief Return the additive identity 0 + 0i.
        /// @return Complex zero.
        static constexpr Complex zero() noexcept { return {}; }

        /// @brief Create a complex number from a real value (helper alias).
        /// @param r  Real part.
        /// @return Complex number r + 0i.
        static constexpr Complex fromReal(FloatType r) noexcept { return { r, 0.0 }; }

		/// @brief Create a complex number from polar coordinates.
		/// @param magnitude  Magnitude (length) of the complex number.
		/// @param angle      Angle (argument) in radians.
		/// @return Complex number with the specified polar representation.
        static constexpr Complex fromPolar(FloatType magnitude, FloatType angle) noexcept
        {
            return { magnitude * ConstexprMath::cos(angle),
                     magnitude * ConstexprMath::sin(angle) };
		}

        /// @brief Return the imaginary unit +i (0 + 1i).
        /// @return Complex representing i.
        static constexpr Complex plus_i() noexcept { return { 0.0, 1.0 }; }

        /// @brief Return the negative imaginary unit -i (0 - 1i).
        /// @return Complex representing -i.
        static constexpr Complex minus_i() noexcept { return { 0.0, -1.0 }; }

        /// @brief Component-wise addition.
        /// @param other  Addend.
        /// @return The sum (this + other).
        constexpr Complex operator+ (Complex other) const noexcept
        {
            // Add real and imaginary parts separately.
            return { re + other.re, im + other.im };
        }

        /// @brief Component-wise subtraction.
        /// @param other  Subtrahend.
        /// @return The difference (this - other).
        constexpr Complex operator- (Complex other) const noexcept
        {
            // Subtract real and imaginary parts separately.
            return { re - other.re, im - other.im };
        }

        /// @brief Unitary minus
        /// @return -this
        constexpr Complex operator-() const noexcept
        {
            return { -re, -im };
        }

        /// @brief Complex multiplication following (a+bi)*(c+di) = (ac - bd) + (ad + bc)i.
        /// @param other  The multiplier.
        /// @return The product.
        constexpr Complex operator* (Complex other) const noexcept
        {
            // Compute using the standard formula to avoid temporaries where possible.
            return { re * other.re - im * other.im, re * other.im + im * other.re };
        }

		/// @brief Scalar multiplication.
		/// @param scalar  The scalar multiplier.
		/// @return The product (this * scalar).
        constexpr Complex operator* (FloatType scalar) const noexcept
        {
            // Scale both parts by the scalar.
            return { re * scalar,  im * scalar };
		}

		/// @brief Scalar division.
		/// @param scalar  The scalar divisor.
		/// @return The quotient (this / scalar).
        constexpr Complex operator/ (FloatType scalar) const noexcept
        {
            // Divide both parts by the scalar.
            return { re / scalar, im / scalar };
		}

        /// @brief Complex division.
        /// @param other The divider (complex).
        /// @return this / other
        constexpr Complex operator/(Complex other) const noexcept
        {
            // Denominator: |other|^2 = c^2 + d^2
            const FloatType Denom =
                other.re * other.re + other.im * other.im;

            // (a+bi)(c-di)
            return
            {
                (re * other.re + im * other.im) / Denom,
                (im * other.re - re * other.im) / Denom
            };
        }

        /// @brief In-place addition (accumulate another complex into this).
        /// @param other  Value to add.
        /// @return A reference to *this after addition (returned by value here for constexpr convenience).
        constexpr Complex operator+= (Complex other) noexcept
        {
            // Mutate both components; returning *this by value is acceptable for constexpr use.
            re += other.re;
            im += other.im;
            return *this;
        }

		/// @brief Compute the complex conjugate.
		/// @return The conjugate of this complex number.
        constexpr Complex conj() const noexcept
        {
            // Conjugate: invert the sign of the imaginary part.
            return { re, -im };
		}

		/// @brief Compute the squared magnitude (norm) of the complex number.
		/// @return The squared norm (re^2 + im^2).
        constexpr FloatType normSquared() const noexcept
        {
            // Return the squared magnitude (norm) of the complex number.
            return re * re + im * im;
		}
    };
}