#pragma once
#include <concepts>
#include "constexpr_trigon.h"

/// @file
/// @brief Constexpr-capable complex number type used throughout the project.
/// @details
/// This header defines `ConstexprMath::Complex<T>`, a small constexpr-friendly
/// complex number type with a minimal set of operations required by the
/// simulator (construction, addition, subtraction, multiplication, in-place
/// accumulation and conjugation).
/// All members and functions are `constexpr` and `noexcept` where appropriate
/// so they can be used in compile-time contexts (e.g., constexpr initialization
/// of gate matrices).

namespace ConstexprMath
{
    /// @brief Simple constexpr complex number for numeric computations.
    /// @tparam T Floating-point type used for real and imaginary parts (e.g. real_t).
    /// @details
    /// The type exposes public data members `re` and `im` and implements lightweight
    /// arithmetic operators in a constexpr-friendly manner. It is intentionally minimal
    /// (no exceptions, no heap allocations).
    template <std::floating_point _FloatType>
    struct Complex
    {
        using FloatType = _FloatType;

        /// @brief Real part of the complex number.
        FloatType re = 0.0;
        /// @brief Imaginary part of the complex number.
        FloatType im = 0.0;

        /// @brief Default constructor producing 0 + 0i.
        constexpr Complex() noexcept = default;

        /// @brief Construct a complex number from a real value (imaginary part = 0).
        /// @param r Real part.
        constexpr Complex(FloatType r) noexcept : re(r), im(0.0) {}

        /// @brief Construct a complex number from real and imaginary parts.
        /// @param r Real component.
        /// @param i Imaginary component.
        constexpr Complex(FloatType r, FloatType i) noexcept : re(r), im(i) {}

        /// @brief Return the additive identity 0 + 0i.
        /// @return Complex zero.
        static constexpr Complex zero() noexcept { return {}; }

        /// @brief Create a complex number from a real value.
        /// @param r Real part.
        /// @return Complex number r + 0i.
        static constexpr Complex fromReal(FloatType r) noexcept { return { r, 0.0 }; }

        /// @brief Create a complex number from polar coordinates.
        /// @param magnitude Magnitude (length) of the number.
        /// @param angle Angle (argument) in radians.
        /// @return Complex value defined by the given polar representation.
        static constexpr Complex fromPolar(FloatType magnitude, FloatType angle) noexcept
        {
            return {
                magnitude * ConstexprMath::cos(angle),
                magnitude * ConstexprMath::sin(angle)
            };
        }

        /// @brief Return +i (0 + 1i).
        /// @return Complex representing i.
        static constexpr Complex plus_i() noexcept { return { 0.0, 1.0 }; }

        /// @brief Return -i (0 - 1i).
        /// @return Complex representing -i.
        static constexpr Complex minus_i() noexcept { return { 0.0, -1.0 }; }

        /// @brief Component-wise addition.
        /// @param other Value to add.
        /// @return The sum (this + other).
        constexpr Complex operator+ (Complex other) const noexcept
        {
            return { re + other.re, im + other.im };
        }

        /// @brief Component-wise subtraction.
        /// @param other Value to subtract.
        /// @return The difference (this - other).
        constexpr Complex operator- (Complex other) const noexcept
        {
            return { re - other.re, im - other.im };
        }

        /// @brief Unary minus.
        /// @return -this.
        constexpr Complex operator-() const noexcept
        {
            return { -re, -im };
        }

        /// @brief Complex multiplication.
        /// @param other Multiplier.
        /// @return The product (this * other).
        constexpr Complex operator* (Complex other) const noexcept
        {
            return { re * other.re - im * other.im,
                     re * other.im + im * other.re };
        }

        /// @brief Scalar multiplication.
        /// @param scalar Scalar multiplier.
        /// @return The product this * scalar.
        constexpr Complex operator* (FloatType scalar) const noexcept
        {
            return { re * scalar, im * scalar };
        }

        /// @brief Scalar division.
        /// @param scalar Scalar divisor.
        /// @return The quotient this / scalar.
        constexpr Complex operator/ (FloatType scalar) const noexcept
        {
            return { re / scalar, im / scalar };
        }

        /// @brief Complex division.
        /// @param other Divisor.
        /// @return this / other.
        constexpr Complex operator/(Complex other) const noexcept
        {
            const FloatType Denom =
                other.re * other.re + other.im * other.im;

            return {
                (re * other.re + im * other.im) / Denom,
                (im * other.re - re * other.im) / Denom
            };
        }

        /// @brief In-place addition.
        /// @param other Addend.
        /// @return *this after mutation (returned by value for constexpr convenience).
        constexpr Complex operator+= (Complex other) noexcept
        {
            re += other.re;
            im += other.im;
            return *this;
        }

        /// @brief Complex conjugate.
        /// @return Conjugate (re, -im).
        constexpr Complex conj() const noexcept
        {
            return { re, -im };
        }

        /// @brief Squared magnitude of the complex number.
        /// @return re² + im².
        constexpr FloatType normSquared() const noexcept
        {
            return re * re + im * im;
        }
    };
}
