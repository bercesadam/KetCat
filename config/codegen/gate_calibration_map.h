// ==========================================================================
//  KetCat AUTOMATICALLY GENERATED OFFLINE CALIBRATION LOOKUP TABLE
//  Execution date: 2026. 06. 28.
// ==========================================================================
#pragma once
#include "core_types.h"
#include "gate_calibration/phase_corrections.h"

#define CALIBRATED_GATES

namespace KetCat
{
    class GateCalibrationTable
    {
    public:
        static constexpr natural_t CalibrationPoints = 9;

        // Exact Newton polynomial interpolator for arbitrary Rx rotations
        static inline PolynomialInterpolator<CalibrationPoints> getRxCalib() noexcept
        {
            constexpr std::array<real_t, CalibrationPoints> x = { 0.0, 0.392699, 0.785398, 1.1781, 1.5708, 1.9635, 2.35619, 2.74889, 3.14159 };
            constexpr std::array<real_t, CalibrationPoints> y = { 0.0, 0.195762, 0.391468, 0.587996, 0.78058, -2.16421, -1.96758, -1.77185, -1.57613 };
            return PolynomialInterpolator<CalibrationPoints>(x, y);
        }

        // Static corrections for the fixed two-qubit Rydberg CPhase (CZ) gate
        static constexpr TwoQubitCalibResult getCPhaseCalib() noexcept
        {
            return { -0.0484911, -0.0484911, 3.15563 };
        }
    };
}


