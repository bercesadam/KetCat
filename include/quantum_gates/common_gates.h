#pragma once
#include "core_types.h"
#include "constexprmath/constexpr_core_functions.h"
#include "constexprmath/constexpr_trigon.h"

/// @file
/// @brief Common single‑ and two‑qubit gate matrices and helper functions.
/// @details
/// Provides commonly used quantum gate matrices (Pauli gates, Hadamard, CNOT,
/// Toffoli, SWAP, rotations) and the `identityMatrix<QBitCount>()` helper for
/// generating identity operators for arbitrary qubit counts.

namespace KetCat::UnitaryGateLib
{
    /// @brief Mathematical constants used for gate definitions.
    constexpr real_t sqrt2     = 1.41421356237309505;
    constexpr real_t inv_sqrt2 = 1.0 / sqrt2;

    /// @brief Produce an identity matrix for `QBitCount` qubits.
    /// @tparam QBitCount Number of qubits.
    /// @return A diagonal identity matrix of size (2^QBitCount × 2^QBitCount).
    template <KetCat::natural_t QBitCount>
    constexpr KetCat::matrix_t<ConstexprMath::pow2(QBitCount)> identityMatrix() noexcept
    {
        constexpr natural_t Dim = ConstexprMath::pow2(QBitCount);

        // Zero-initialize and set diagonal elements to 1.
        matrix_t<Dim> Identity = {};
        for (natural_t i = 0; i < Dim; ++i)
        {
            Identity[i][i] = complex_t::fromReal(1.0);
        }
        return Identity;
    }

    /// @brief Single-qubit identity matrix.
    constexpr KetCat::matrix_t<2> I = identityMatrix<1U>();

    /// @brief Hadamard gate:
    ///        H = (1/sqrt(2)) * [[1, 1], [1, −1]]
    constexpr KetCat::matrix_t<2> H = { {
        { complex_t(1.0 / sqrt2, 0.0),  complex_t(1.0 / sqrt2, 0.0) },
        { complex_t(1.0 / sqrt2, 0.0),  complex_t(-1.0 / sqrt2, 0.0) }
    } };

    /// @brief Pauli‑X gate (NOT).
    constexpr KetCat::matrix_t<2> X = { {
        { complex_t(0.0, 0.0), complex_t(1.0, 0.0) },
        { complex_t(1.0, 0.0), complex_t(0.0, 0.0) }
    } };

    /// @brief Pauli‑Y gate.
    constexpr KetCat::matrix_t<2> Y = { {
        { complex_t(0.0, 0.0),  complex_t(0.0, -1.0) },
        { complex_t(0.0, 1.0),  complex_t(0.0, 0.0) }
    } };

    /// @brief Pauli‑Z gate.
    constexpr KetCat::matrix_t<2> Z = { {
        { complex_t(1.0, 0.0),  complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0),  complex_t(-1.0, 0.0) }
    } };

    /// @brief Two‑qubit SWAP gate.
    constexpr KetCat::matrix_t<4> SWAP = { {
        { complex_t(1.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(1.0, 0.0), complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(1.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(1.0, 0.0) }
    } };

    /// @brief CNOT gate (control is the first qubit).
    constexpr KetCat::matrix_t<4> CX = { {
        { complex_t(1.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(1.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(0.0, 0.0), complex_t(1.0, 0.0), complex_t(0.0, 0.0) },
        { complex_t(0.0, 0.0), complex_t(1.0, 0.0), complex_t(0.0, 0.0), complex_t(0.0, 0.0) }
    } };
    constexpr auto CNOT = CX;

    /// @brief Toffoli gate (CCNOT / CCX), control-first.
    constexpr KetCat::matrix_t<8> CCX = { {
        { complex_t(1.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(1.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(1.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(1.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(1.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(1.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(1.0,0.0) },

        { complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(0.0,0.0),
          complex_t(0.0,0.0), complex_t(0.0,0.0), complex_t(1.0,0.0), complex_t(0.0,0.0) }
    } };
    constexpr auto TOFFOLI = CCX;

    /// @brief Rotation around the Y‑axis by angle θ.
    /// @details Implements R_y(θ) = cos(θ/2) I - i sin(θ/2) Y
    /// @param theta Rotation angle in radians.
    /// @return 2×2 rotation matrix.
    constexpr KetCat::matrix_t<2> RotationY(real_t theta) noexcept
    {
        const real_t halfTheta = theta / 2.0;
        return KetCat::matrix_t<2>{ {
            { complex_t(ConstexprMath::cos(halfTheta), 0.0),
              complex_t(-ConstexprMath::sin(halfTheta), 0.0) },
            { complex_t(ConstexprMath::sin(halfTheta), 0.0),
              complex_t(ConstexprMath::cos(halfTheta), 0.0) }
        } };
    }

} // namespace KetCat::QCC::Gates
