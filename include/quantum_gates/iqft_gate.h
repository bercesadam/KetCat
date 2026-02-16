#pragma once
#include "common_gates.h"
#include "constexprmath/constexpr_trigon.h"

namespace KetCat::QCC::Gates
{
    /// @brief Construct the inverse Quantum Fourier Transform (QFT†) dense matrix for QBitCount qubits.
    /// @tparam QBitCount Number of qubits (matrix dimension = 2^QBitCount).
    /// @return A unitary matrix implementing the inverse QFT:
    ///         U[j,k] = (1/√N) · exp(−2π i j k / N), where N = 2^QBitCount.
    /// @note Constexpr‑friendly and uses the project’s constexpr trigonometric helpers
    ///       for compile‑time matrix construction.
    template<dimension_t QBitCount>
    constexpr matrix_t<ConstexprMath::pow2(QBitCount)>
        make_IQFT_matrix() noexcept
    {
        constexpr dimension_t Dim = ConstexprMath::pow2(QBitCount);
        constexpr real_t InvSqrtDim = 1.0 / ConstexprMath::sqrt(static_cast<real_t>(Dim));

        matrix_t<Dim> IQFTMatrix{};

        for (dimension_t j = 0; j < Dim; ++j)
        {
            for (dimension_t k = 0; k < Dim; ++k)
            {
                // angle = 2π · j · k / N
                const real_t Angle =
                    2.0 * ConstexprMath::Pi * static_cast<real_t>(j * k)
                    / static_cast<real_t>(Dim);

                // Inverse QFT uses negative exponent:
                // exp(−i · angle) = cos(angle) − i · sin(angle)
                const real_t re = ConstexprMath::cos(Angle);
                const real_t im = -ConstexprMath::sin(Angle);

                IQFTMatrix[j][k] = cplx_t(
                    InvSqrtDim * re,
                    InvSqrtDim * im
                );
            }
        }

        return IQFTMatrix;
    }
}
