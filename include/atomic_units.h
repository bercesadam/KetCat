#pragma once
#include "core_types.h"


namespace KetCat::Units
{
    /// -------------------------------------------------------------------------
    /// Atomic units  according to CODATA  
    /// -------------------------------------------------------------------------

    /// Atomic unit of time (ℏ / Eh)
    constexpr real_t AtomicTimeToSeconds = 2.4188843265857e-17; // s

    /// Atomic unit of intensity (corresponding to 1 a.u. electric field squared)
    constexpr real_t AtomicIntensityToWcm2 = 3.509445219e16; // W/cm²

    /// Speed of light in vacuum
    constexpr real_t SpeedOfLightMps = 299792458.0; // m/s

    // 1 a.u. angular frequency in Hz: omega_au / AtomicTimeToSeconds
    constexpr real_t AuOmegaToHz = 1.0 / AtomicTimeToSeconds;

    
    // -------------------------------------------------------------------------
    // Angular frequency conversions
    // -------------------------------------------------------------------------

    /// Conversion factor: ω[a.u.] → ω[rad/s]
    constexpr real_t AtomicOmegaToRadPerSecond = 1.0 / AtomicTimeToSeconds;

    /// Convert angular frequency from rad/s to atomic units
    constexpr real_t omegaAuFromRadPerSecond(real_t omegaRadPerSecond) noexcept
    {
        return omegaRadPerSecond * AtomicTimeToSeconds;
    }

    /// Convert angular frequency from atomic units to rad/s
    constexpr real_t radPerSecondFromOmegaAu(real_t omegaAu) noexcept
    {
        return omegaAu / AtomicTimeToSeconds;
    }

    /// Convert frequency (cycles/s, Hz) to angular frequency in atomic units
    /// ω = 2π f
    constexpr real_t omegaAuFromHz(real_t frequencyHz) noexcept
    {
        return ConstexprMath::Pi * 2 * frequencyHz * AtomicTimeToSeconds;
    }

    /// Convert angular frequency in atomic units to frequency in Hz (cycles/s)
    constexpr real_t hzFromOmegaAu(real_t omegaAu) noexcept
    {
        const real_t omegaSI = omegaAu / AtomicTimeToSeconds; // rad/s
        return omegaSI / ConstexprMath::Pi * 2;                              // Hz
    }


    // -------------------------------------------------------------------------
    // Wavelength ↔ angular frequency
    // -------------------------------------------------------------------------

    /// Convert angular frequency (a.u.) to wavelength (nm)
    /// λ = 2π c / ω
    constexpr real_t wavelengthNmFromOmegaAu(real_t omegaAu) noexcept
    {
        const real_t omegaSI = omegaAu / AtomicTimeToSeconds; // rad/s
        const real_t wavelengthM =
            (ConstexprMath::Pi * 2 * SpeedOfLightMps) / omegaSI;

        return wavelengthM * 1e9; // nm
    }

    /// Convert wavelength (nm) to angular frequency (a.u.)
    /// ω = 2π c / λ
    constexpr real_t omegaAuFromWavelengthNm(real_t wavelengthNm) noexcept
    {
        const real_t wavelengthM = wavelengthNm * 1e-9;
        const real_t omegaSI =
            (ConstexprMath::Pi * 2 * SpeedOfLightMps) / wavelengthM;

        return omegaSI * AtomicTimeToSeconds;
    }


    // -------------------------------------------------------------------------
    // Electric field ↔ intensity
    // -------------------------------------------------------------------------

    /// Convert electric field amplitude (a.u.) to intensity (W/cm²)
    /// I = I_au · E²
    constexpr real_t intensityWcm2FromFieldAu(real_t fieldAmplitudeAu) noexcept
    {
        return AtomicIntensityToWcm2
             * fieldAmplitudeAu
             * fieldAmplitudeAu;
    }

    /// Convert intensity (W/cm²) to electric field amplitude (a.u.)
    constexpr real_t fieldAuFromIntensityWcm2(real_t intensityWcm2) noexcept
    {
        return ConstexprMath::sqrt(
            intensityWcm2 / AtomicIntensityToWcm2
        );
    }
}
