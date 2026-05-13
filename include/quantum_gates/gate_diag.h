#pragma once
#include "hilbert_space/matrix.h"
#include "hilbert_space/state_vector.h"
#include "quantum_gates/gates.h"


namespace KetCat
{
    /// @brief Diagnostics and characterization utilities for quantum gates.
    ///
    /// @details
    ///   This utility extracts:
    ///
    ///     • Effective gate matrices from propagated basis states
    ///     • Global phase offsets
    ///     • Average gate fidelity
    ///     • Process overlap metrics
    ///
    ///   The extracted effective gate is:
    ///
    ///       U_eff = [ |ψ₀_out⟩  |ψ₁_out⟩ ]
    ///
    ///   where:
    ///
    ///       |ψ₀_out⟩ = U |0⟩
    ///       |ψ₁_out⟩ = U |1⟩
    ///
    ///   The ideal comparison gate may be any unitary matrix.
    template<natural_t Dim>
    class GateDiagnostic
    {
    public:

        using Matrix = Matrix<Dim>;

        /// @brief Construct effective gate matrix from output basis states.
        ///
        /// @param basisOutputs
        ///   Array containing propagated computational basis outputs.
        ///
        /// @return
        ///   Effective unitary approximation.
        ///
        /// @details
        ///   Example for single-qubit systems:
        ///
        ///       basisOutputs[0] = U|0⟩
        ///       basisOutputs[1] = U|1⟩
        ///
        ///   Then:
        ///
        ///            [ ψ₀₀ ψ₁₀ ]
        ///       U =  [ ψ₀₁ ψ₁₁ ]
        ///
        template<typename HilbertSpace>
        static constexpr Matrix buildEffectiveGate(
            const std::array<StateVector<HilbertSpace>, Dim>& basisOutputs) noexcept
        {
            Matrix Result;
            Result.setZero();

            for (natural_t col = 0; col < Dim; ++col)
            {
                for (natural_t row = 0; row < Dim; ++row)
                {
                    Result(row, col) =
                        basisOutputs[col][row];
                }
            }

            return Result;
        }

        /// @brief Compute global phase difference between two gates.
        ///
        /// @param effective
        ///   Numerically extracted gate.
        ///
        /// @param ideal
        ///   Target ideal gate.
        ///
        /// @return
        ///   Estimated global phase in radians.
        ///
        /// @details
        ///   Uses:
        ///
        ///       φ = arg(Tr(U_ideal† U_eff))
        ///
        static constexpr real_t globalPhase(
            const Matrix& effective,
            const Matrix& ideal) noexcept
        {
            const Matrix Overlap =
                ideal.dagger() * effective;

            const complex_t Tr =
                Overlap.trace();

            return ConstexprMath::atan2(
                Tr.im,
                Tr.re);
        }

        /// @brief Remove global phase from a gate matrix.
        ///
        /// @param gate
        ///   Input gate matrix.
        ///
        /// @param phase
        ///   Phase angle in radians.
        ///
        /// @return
        ///   Phase-corrected matrix.
        static constexpr Matrix removeGlobalPhase(
            const Matrix& gate,
            real_t phase) noexcept
        {
            const complex_t PhaseFactor =
            {
                ConstexprMath::cos(-phase),
                ConstexprMath::sin(-phase)
            };

            return gate * PhaseFactor;
        }

        /// @brief Compute process overlap between two gates.
        ///
        /// @details
        ///   Computes:
        ///
        ///       Tr(U_ideal† U_eff)
        ///
        static constexpr complex_t processOverlap(
            const Matrix& effective,
            const Matrix& ideal) noexcept
        {
            return (ideal.dagger() * effective).trace();
        }

        /// @brief Compute average gate fidelity.
        ///
        /// @details
        ///   Uses the standard unitary fidelity:
        ///
        ///                    |Tr(U†V)|² + d
        ///       F_avg = -------------------------
        ///                   d(d + 1)
        ///
        ///   where:
        ///
        ///       d = Hilbert space dimension
        ///
        static constexpr real_t averageGateFidelity(
            const Matrix& effective,
            const Matrix& ideal) noexcept
        {
            const complex_t Overlap =
                processOverlap(effective, ideal);

            const real_t Abs2 =
                Overlap.normSquared();

            constexpr real_t d =
                static_cast<real_t>(Dim);

            return (Abs2 + d) / (d * (d + 1.0));
        }

        /// @brief Compute operator norm deviation from target gate.
        ///
        /// @details
        ///   Computes Frobenius norm:
        ///
        ///       ||U_eff - U_ideal||
        ///
        static constexpr real_t gateErrorNorm(
            const Matrix& effective,
            const Matrix& ideal) noexcept
        {
            real_t Sum = 0.0;

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    const complex_t Diff =
                        effective(i, j) - ideal(i, j);

                    Sum += Diff.normSquared();
                }
            }

            return ConstexprMath::sqrt(Sum);
        }

        /// @brief Check approximate unitarity.
        ///
        /// @param epsilon
        ///   Allowed numerical tolerance.
        static constexpr bool isApproximatelyUnitary(
            const Matrix& gate,
            real_t epsilon = 1E-6) noexcept
        {
            const Matrix Identity =
                Matrix::identity();

            const Matrix Product =
                gate * gate.dagger();

            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    const complex_t Diff =
                        Product(i, j) - Identity(i, j);

                    if (Diff.normSquared() >
                        epsilon * epsilon)
                    {
                        return false;
                    }
                }
            }

            return true;
        }
    };
}