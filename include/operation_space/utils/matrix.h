#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Fixed-size complex square matrix with constexpr linear algebra utilities.
    ///
    /// @tparam Dim
    ///   Matrix dimension.
    ///
    /// @details
    ///   Provides lightweight constexpr-compatible matrix operations
    ///   used throughout the quantum simulation framework:
    ///
    ///     • Matrix multiplication
    ///     • Matrix-vector multiplication
    ///     • Hermitian conjugation (dagger)
    ///     • Trace
    ///     • Identity construction
    ///     • Transpose
    ///
    ///   Intended usage:
    ///
    ///     • Quantum gates
    ///     • Effective logical operators
    ///     • Density matrices
    ///     • Tensor-product extracted subspaces
    ///
    ///   Storage layout:
    ///
    ///       m[row][column]
    template<natural_t Dim>
    struct Matrix
    {
        /// @brief Underlying matrix storage.
        matrix_t<Dim> m{};

        /// INDEXING /////////////////////////////////////////////////////////

        constexpr std::array<complex_t, Dim>&
            operator[](natural_t row) noexcept
        {
            return m[row];
        }

        constexpr const std::array<complex_t, Dim>&
            operator[](natural_t row) const noexcept
        {
            return m[row];
        }

        /// BASIC UTILITIES //////////////////////////////////////////////////

        /// @brief Set all elements to zero.
        constexpr void setZero() noexcept
        {
            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    m[i][j] = complex_t::zero();
                }
            }
        }

        /// @brief Create identity matrix.
        static constexpr Matrix identity() noexcept
        {
            Matrix Result;
            Result.setZero();

            for (natural_t i = 0; i < Dim; ++i)
            {
                Result[i][i] =
                    complex_t::fromReal(1.0);
            }

            return Result;
        }

        /// @brief Compute matrix trace.
        ///
        /// @details
        ///     Tr(A) = Σ Aᵢᵢ
        constexpr complex_t trace() const noexcept
        {
            complex_t Result =
                complex_t::zero();

            for (natural_t i = 0; i < Dim; ++i)
            {
                Result += m[i][i];
            }

            return Result;
        }

        /// @brief Compute transpose.
        ///
        /// @details
        ///     Aᵀ[i,j] = A[j,i]
        constexpr Matrix transpose() const noexcept
        {
            Matrix Result;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result[i][j] = m[j][i];
                }
            }

            return Result;
        }

        /// @brief Compute Hermitian conjugate.
        ///
        /// @details
        ///     A† = (Aᵀ)*
        constexpr Matrix dagger() const noexcept
        {
            Matrix Result;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result[i][j] =
                        m[j][i].conj();
                }
            }

            return Result;
        }

        /// MATRIX ALGEBRA ///////////////////////////////////////////////////

        /// @brief Matrix multiplication.
        ///
        /// @details
        ///     C = A · B
        constexpr Matrix
            operator*(const Matrix& rhs) const noexcept
        {
            Matrix Result;
            Result.setZero();

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    complex_t Sum =
                        complex_t::zero();

                    for (natural_t k = 0; k < Dim; ++k)
                    {
                        Sum +=
                            m[i][k]
                            * rhs[k][j];
                    }

                    Result[i][j] = Sum;
                }
            }

            return Result;
        }

        /// @brief Matrix addition.
        constexpr Matrix
            operator+(const Matrix& rhs) const noexcept
        {
            Matrix Result;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result[i][j] =
                        m[i][j]
                        + rhs[i][j];
                }
            }

            return Result;
        }

        /// @brief Matrix subtraction.
        constexpr Matrix
            operator-(const Matrix& rhs) const noexcept
        {
            Matrix Result;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result[i][j] =
                        m[i][j]
                        - rhs[i][j];
                }
            }

            return Result;
        }

        /// @brief Scalar multiplication.
        constexpr Matrix
            operator*(complex_t scalar) const noexcept
        {
            Matrix Result;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result[i][j] =
                        m[i][j] * scalar;
                }
            }

            return Result;
        }

        /// @brief Matrix-vector multiplication.
        ///
        /// @details
        ///     |ψ'⟩ = A |ψ⟩
        constexpr state_vector_t<Dim>
            operator*(const state_vector_t<Dim>& vec) const noexcept
        {
            state_vector_t<Dim> Result{};

            for (natural_t i = 0; i < Dim; ++i)
            {
                complex_t Sum =
                    complex_t::zero();

                for (natural_t j = 0; j < Dim; ++j)
                {
                    Sum +=
                        m[i][j] * vec[j];
                }

                Result[i] = Sum;
            }

            return Result;
        }

        /// QUANTUM HELPERS /////////////////////////////////////////////////

        /// @brief Check approximate unitarity.
        ///
        /// @param tolerance
        ///   Numerical tolerance.
        ///
        /// @details
        ///   Verifies:
        ///
        ///       U†U ≈ I
        constexpr bool
            isUnitary(real_t tolerance = 1E-9) const noexcept
        {
            const Matrix Product =
                dagger() * (*this);

            const Matrix Identity =
                Matrix::identity();

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    const complex_t Diff =
                        Product[i][j]
                        - Identity[i][j];

                    if (Diff.normSquared()
                > tolerance * tolerance)
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        /// @brief Frobenius norm squared.
        ///
        /// @details
        ///     ||A||² = Σ |Aᵢⱼ|²
        constexpr real_t
            frobeniusNormSquared() const noexcept
        {
            real_t Result = 0.0;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Result +=
                        m[i][j].normSquared();
                }
            }

            return Result;
        }
    };
}