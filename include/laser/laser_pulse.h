#pragma once
#include "core_types.h"
#include "hamiltonian/rabi_drive_hamiltonian.h"


namespace KetCat
{
    /// @brief Converts laser parameters from SI units to Atomic Units (a.u.).
    ///
    /// @details
    /// This utility provides conversion constants and logic to bridge the gap between 
    /// laboratory units (nm, W/cm²) and the internal atomic units used by the solvers.
    ///
    /// Constants used (CODATA):
    /// - τ₀ ≈ 2.4189e-17 s (Atomic unit of time)
    /// - I₀ ≈ 3.5094e16 W/cm² (Atomic unit of intensity)
    struct SiLaserPulse
    {
        // Conversion constants
        static constexpr real_t AU_TIME = 2.4188843265857e-17;       ///< [s]
        static constexpr real_t AU_INTENSITY = 3.509445219e16;       ///< [W/cm²]
        static constexpr real_t SPEED_OF_LIGHT_SI = 299792458.0;     ///< [m/s]

        real_t omega_au;      ///< Angular frequency [Eₕ/ħ]
        real_t amplitude_au;  ///< Peak electric field amplitude [Eₕ/(e·a₀)]

        /// @brief Calculates atomic parameters from wavelength and intensity.
        ///
        /// @param wavelength_nm  Laser wavelength λ [nm]
        /// @param intensity_Wcm2 Laser peak intensity I [W/cm²]
        ///
        /// @details
        /// 1. Frequency: ω_si = 2πc / λ
        /// 2. Field:     ε₀ = √(I / I_au)
        static /*constexpr*/ SiLaserPulse fromSi(real_t wavelength_nm, real_t intensity_Wcm2) noexcept
        {
            SiLaserPulse pulse{};
            pulse.omega_au = KetCat::Units::omegaAuFromWavelengthNm(wavelength_nm);
            pulse.amplitude_au = KetCat::Units::fieldAuFromIntensityWcm2(intensity_Wcm2);

			std::cout << "Converted SI parameters to AU: " << std::endl;
			std::cout << "  Wavelength (nm): " << wavelength_nm << " -> Angular Frequency (a.u.): " << pulse.omega_au << std::endl;
			std::cout << "  Intensity (W/cm²): " << intensity_Wcm2 << " -> Field Amplitude (a.u.): " << pulse.amplitude_au << std::endl;

            return pulse;
        }
    };

    /// @brief Factory for instantiating Hamiltonians using SI laser parameters.
    ///
    /// @tparam LevelCount The number of levels in the energy space.
    template<natural_t LevelCount>
    class LaserHamiltonianBuilder
    {
    public:
        /// @brief Builds a RwaRabiHamiltonian by converting SI inputs to AU.
        ///
        /// @param energies       Energies in AU [Eₕ].
        /// @param dipoleMatrix   Radial dipole moments in AU [e·a₀].
        /// @param wavelength_nm  Input wavelength [nm].
        /// @param intensity_Wcm2 Input intensity [W/cm²].
        /// @param referenceLevel Index defining the rotating frame (default = 0).
        static constexpr auto build(
            const std::array<real_t, LevelCount>& energies,
            const matrix_t<LevelCount>& dipoleMatrix,
            const real_t wavelength_nm,
            const real_t intensity_Wcm2,
            const natural_t referenceLevel = 0) noexcept
        {
            // Perform SI to AU conversion
            const SiLaserPulse pulse = SiLaserPulse::fromSi(wavelength_nm, intensity_Wcm2);

            // Forward to the solver-compatible Hamiltonian class
            return RwaRabiHamiltonian<LevelCount>(
                energies,
                dipoleMatrix,
                pulse.omega_au,
                pulse.amplitude_au,
                referenceLevel
            );
        }
    };
}