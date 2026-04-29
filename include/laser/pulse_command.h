#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Quantization axes for Bloch sphere rotations.
    enum class RotationAxis
    {
        X, ///< Rotation around the X-axis (in-phase drive).
        Y, ///< Rotation around the Y-axis (quadrature drive, 90-degree shift).
        Z  ///< Rotation around the Z-axis (longitudinal phase shift/detuning).
    };

    /// @brief High-level command for a quantum gate operation.
    struct PulseCommand
    {
        RotationAxis m_axis;        ///< Target axis of rotation.
        real_t m_rotationAngleRad;  ///< Rotation angle θ in radians.
    };

}