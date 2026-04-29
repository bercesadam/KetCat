#pragma once
#include <tuple>
#include "atomic_units.h"
#include "laser/laser.h"
#include "systems/neutral_atom_manifold.h"


namespace KetCat
{
    /// @brief Time step for Time-Dependent Schrödinger Equation (TDSE) integration in atomic units (a.u.).
    static constexpr real_t TimeStepAu = 50;

    struct STIRAPConfig
    {
        real_t m_Level1Energy;
        real_t m_Level2Energy;
        real_t m_Level3Energy;

        real_t m_Mu12;
        real_t m_Mu23;

        real_t m_rabiFrequency_Hz;
        real_t m_pumpPhase;
        real_t m_stokesDetuning;
        real_t m_targetTheta;
    };


    template <NeutralAtomTypeConfig Config>
    class STIRAPLaser
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;

		STIRAPConfig m_config;     //< Externally provided configuration for the STIRAP pulse sequence, including target rotation angle and detuning.

		real_t m_Sigma;			   ///< Width of the Gaussian pulses (in a.u.), related to the pulse duration.
		real_t m_Omega0;           ///< Peak Rabi frequency (in a.u.) for the pulse sequence.

		real_t m_StokesTimeCenter; ///< Center time for the Stokes pulse (in a.u.).
		real_t m_PumpTimeCenter;   ///< Center time for the Pump pulse (in a.u.).
		real_t m_TransitionTimeLimit; ///< Time limit for truncating the pulse sequence in fractional STIRAP (in a.u.).

		real_t m_StokesOmega;       ///< Rabi frequency for the Stokes pulse (in a.u.).
		real_t m_PumpOmega;         ///< Rabi frequency for the Pump pulse (in a.u.).

        /// @brief Simple Gaussian temporal profile for the Rabi frequency.
        ///
        /// @param t      Current time.
        /// @param t0     Peak position (mean).
        /// @param sigma  Width of the pulse (standard deviation).
        /// @return Pulse amplitude at time 't' relative to the peak.
        real_t gaussian(const real_t t, const real_t t0, const real_t sigma) const
        {
            real_t X = (t - t0) / sigma;
            return std::exp(-X * X);
        }

    public:
        constexpr STIRAPLaser(const STIRAPConfig& config)
            : m_config(config)
        {
            // Target Rabi frequency parameters for defining the pulse area.
            m_Omega0 = Units::omegaAuFromHz(config.m_rabiFrequency_Hz);
            m_Sigma = 50.0 / m_Omega0;

            // Center-times for the delayed pulse sequence (Stokes precedes Pump).
            m_StokesTimeCenter = 4.0 * m_Sigma;
            m_PumpTimeCenter = 6.0 * m_Sigma;
            real_t DeltaTime = m_PumpTimeCenter - m_StokesTimeCenter;

			
			real_t Ratio = 1.0; // Default value for θ = π (full population transfer)

            // Calculate Fractional STIRAP limit for the desired rotation angle θ.
			// Protecting against edge cases where θ is very close to 0 or π to avoid singulariy for the tan function
            if (ConstexprMath::abs(config.m_targetTheta - ConstexprMath::Pi) > 1E-9)
            {
                Ratio = ConstexprMath::tan(config.m_targetTheta / 2.0);
            }

			// For zero rotation (θ = 0), skip the pulse sequence entirely
            if (config.m_targetTheta <= 0.0)
            {
                m_TransitionTimeLimit = 0.0;
            }
			// Full population transfer (θ = π)
            else if (config.m_targetTheta >= ConstexprMath::Pi)
            {
                m_TransitionTimeLimit = 10.0 * m_Sigma;
            }
			// Fractional STIRAP (0 < θ < π)
			// t_limit = ((t_stokes + t_pump) / 2) + (σ² / Δt) * ln(tan(θ/2))
            else
            {
                m_TransitionTimeLimit = ((m_PumpTimeCenter + m_StokesTimeCenter) / 2.0) +
                    ((m_Sigma * m_Sigma) / DeltaTime) * ConstexprMath::log(Ratio);
            }

            // Extract relevant parameters from the manifold for later use in the time dependent Laser generation.
            m_PumpOmega = config.m_Level2Energy - config.m_Level1Energy; // Resonant with the first transition
			m_StokesOmega = (config.m_Level3Energy - config.m_Level2Energy) + config.m_stokesDetuning; // Resonant with the second transition + optional Z-detuning
        }

        std::tuple<LaserPulse, LaserPulse> operator()(real_t time) const
        {
            real_t OmegaStokes = m_Omega0 * gaussian(time, m_StokesTimeCenter, m_Sigma);
            real_t OmegaPump = m_Omega0 * gaussian(time, m_PumpTimeCenter, m_Sigma);

            // Define the Rotating Wave Approximation (RWA) lasers.
            LaserPulse Stokes{ };
            LaserPulse Pump{ };

            // Laser 0 (Pump): Resonant with the first transition, phase-controlled.
            Pump = { m_PumpOmega, OmegaPump / m_config.m_Mu12, m_config.m_pumpPhase };

            // Laser 1 (Stokes): Resonant with the second transition + optional Z-detuning.
            Stokes = { m_StokesOmega, OmegaStokes / m_config.m_Mu23, 0.0 };

            return { Pump, Stokes };
        }

        real_t getTransitionTimeLimit() const
        {
            return m_TransitionTimeLimit;
		}
    };
}
