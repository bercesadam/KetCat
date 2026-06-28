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
            constexpr std::array<real_t, CalibrationPoints> x = {0.0, .392699, 0.785398, 1.1781, 1.5708, 1.9635, 2.35619, 2.74889, 3.14159};
            constexpr std::array<real_t, CalibrationPoints> y = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            return PolynomialInterpolator<CalibrationPoints>(x, y);
        }

        // Static corrections for the fixed two-qubit Rydberg CPhase (CZ) gate
        static constexpr TwoQubitCalibResult getCPhaseCalib() noexcept
        {
            return { 1.0968, 1.0968, -3.13655 };
        }
    };
}
