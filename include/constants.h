#pragma once
#include "core_types.h"


namespace KetCat::Constants
{
    constexpr real_t Hbar = 1.0; // Planck's constant divided by 2π, in atomic units
    constexpr real_t ElectronMass = 1.0; // Electron mass, in atomic units
    constexpr real_t ElementaryCharge = 1.0; // Elementary charge, in atomic units
    constexpr real_t BohrRadius = 1.0; // Bohr radius, in atomic units
}

namespace KetCat::Units
{
    // Atomic units (CODATA)    
    constexpr real_t TwoPi                      = ConstexprMath::Pi * 2;

    /// Atomic unit of time (ℏ / Eh)
    inline constexpr real_t atomicTimeToSeconds = 2.4188843265857e-17; // s

    /// Atomic unit of intensity (corresponding to 1 a.u. electric field squared)
    inline constexpr real_t atomicIntensityToWcm2 = 3.509445219e16; // W/cm²

    /// Speed of light in vacuum
    inline constexpr real_t speedOfLightMps = 299792458.0; // m/s

    // 1 a.u. angular frequency in Hz: omega_au / atomicTimeToSeconds
    inline constexpr real_t auOmega_to_Hz = 1.0 / atomicTimeToSeconds;

    // -------------------------------------------------------------------------
    // Angular frequency conversions
    // -------------------------------------------------------------------------

    /// Conversion factor: ω[a.u.] → ω[rad/s]
    inline constexpr real_t atomicOmegaToRadPerSecond =
        1.0 / atomicTimeToSeconds;

    /// Convert angular frequency from rad/s to atomic units
    inline real_t omegaAuFromRadPerSecond(real_t omegaRadPerSecond) noexcept
    {
        return omegaRadPerSecond * atomicTimeToSeconds;
    }

    /// Convert angular frequency from atomic units to rad/s
    inline real_t radPerSecondFromOmegaAu(real_t omegaAu) noexcept
    {
        return omegaAu / atomicTimeToSeconds;
    }

    /// Convert frequency (cycles/s, Hz) to angular frequency in atomic units
    /// ω = 2π f
    inline real_t omegaAuFromHz(real_t frequencyHz) noexcept
    {
        return TwoPi * frequencyHz * atomicTimeToSeconds;
    }

    /// Convert angular frequency in atomic units to frequency in Hz (cycles/s)
    inline real_t hzFromOmegaAu(real_t omegaAu) noexcept
    {
        const real_t omegaSI = omegaAu / atomicTimeToSeconds; // rad/s
        return omegaSI / TwoPi;                              // Hz
    }


    // -------------------------------------------------------------------------
    // Wavelength ↔ angular frequency
    // -------------------------------------------------------------------------

    /// Convert angular frequency (a.u.) to wavelength (nm)
    /// λ = 2π c / ω
    inline real_t wavelengthNmFromOmegaAu(real_t omegaAu) noexcept
    {
        const real_t omegaSI = omegaAu / atomicTimeToSeconds; // rad/s
        const real_t wavelengthM =
            (TwoPi * speedOfLightMps) / omegaSI;

        return wavelengthM * 1e9; // nm
    }

    /// Convert wavelength (nm) to angular frequency (a.u.)
    /// ω = 2π c / λ
    inline real_t omegaAuFromWavelengthNm(real_t wavelengthNm) noexcept
    {
        const real_t wavelengthM = wavelengthNm * 1e-9;
        const real_t omegaSI =
            (TwoPi * speedOfLightMps) / wavelengthM;

        return omegaSI * atomicTimeToSeconds;
    }


    // -------------------------------------------------------------------------
    // Electric field ↔ intensity
    // -------------------------------------------------------------------------

    /// Convert electric field amplitude (a.u.) to intensity (W/cm²)
    /// I = I_au · E²
    inline real_t intensityWcm2FromFieldAu(real_t fieldAmplitudeAu) noexcept
    {
        return atomicIntensityToWcm2
             * fieldAmplitudeAu
             * fieldAmplitudeAu;
    }

    /// Convert intensity (W/cm²) to electric field amplitude (a.u.)
    inline real_t fieldAuFromIntensityWcm2(real_t intensityWcm2) noexcept
    {
        return ConstexprMath::sqrt(
            intensityWcm2 / atomicIntensityToWcm2
        );
    }
}
