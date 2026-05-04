#pragma once
#include <tuple>
#include "atomic_units.h"
#include "laser/laser.h"
#include "systems/neutral_atom_manifold.h"


namespace KetCat
{
    /**
        * @brief Control protocols for 3-level ladder systems.
        */
    enum class TwoPhotonPulseProtocol
    {
        STIRAP,         ///< Stokes then Pump: Robust transfer 0 -> 2
        InvertedSTIRAP, ///< Pump then Stokes: Robust transfer 2 -> 0
        Simultaneous    ///< Raman-resonant: Precise rotation by target theta
    };

    /**
        * @brief Physical configuration for two-photon transition.
        */
    struct TwoPhotonConfig
    {
        real_t m_Level1Energy;      ///< E0 (Ground state)
        real_t m_Level2Energy;      ///< E1 (Intermediate state, e.g., 7p)
        real_t m_Level3Energy;      ///< E2 (Target state, e.g., 7s)
        real_t m_Mu12;              ///< Transition dipole <0|μ|1>
        real_t m_Mu23;              ///< Transition dipole <1|μ|2>

        real_t m_peakRabiFrequency; ///< Peak Ω_0 (in a.u.) for both pulses
        real_t m_pumpPhase;         ///< Phase φ relative to Stokes
        real_t m_commonDetuning;    ///< Δ: Detuning from level 2
        real_t m_targetTheta;       ///< Rotation angle on the 0-2 Bloch sphere (Pi = full transfer)

        TwoPhotonPulseProtocol m_protocol = TwoPhotonPulseProtocol::STIRAP;
    };

    /**
        * @brief Precise Laser Pulse Generator for Ladder Systems.
        * @details Provides time-dependent LaserPulse objects for TDSE solvers.
        */
    class GenericTwoPhotonLaser
    {
        TwoPhotonConfig m_config;
        real_t m_Sigma;
        real_t m_tP, m_tS;      ///< Centers of Pump and Stokes pulses
        real_t m_TimeLimit;     ///< Total duration for the solver to run
        real_t m_omegaP;        ///< Physical frequency for Pump laser
        real_t m_omegaS;        ///< Physical frequency for Stokes laser

        static constexpr real_t gaussian(real_t t, real_t t0, real_t sigma) noexcept
        {
            if (sigma < 1e-12) return 0.0;
            const real_t x = (t - t0) / sigma;
            return std::exp(-x * x);
        }

    public:
        constexpr GenericTwoPhotonLaser(const TwoPhotonConfig& config)
            : m_config(config)
        {
            // Initial frequencies for a ladder: E0 -> E1 -> E2
            // To maintain two-photon resonance: (E2 - E0) = hBar * (omegaP + omegaS)
            // If Pump is detuned by +Delta, Stokes must be detuned by -Delta.
            m_omegaP = (m_config.m_Level2Energy - m_config.m_Level1Energy) + m_config.m_commonDetuning;
            m_omegaS = (m_config.m_Level3Energy - m_config.m_Level2Energy) - m_config.m_commonDetuning;

            if (m_config.m_protocol == TwoPhotonPulseProtocol::Simultaneous)
            {
                setupRamanParameters();
            }
            else
            {
                setupStirapParameters();
            }
        }

    private:
        void setupRamanParameters()
        {
            // Raman effective Rabi frequency: Ω_eff(t) = (Ωp(t) * Ωs(t)) / (2 * Δ)
            // For simultaneous identical Gaussians: Ω_eff(t) = peakEff * exp(-2 * (t-t0)^2 / σ^2)
            real_t peakOmegaEff;
            if (std::abs(m_config.m_commonDetuning) < 1e-10)
            {
                // Resonant case (not strictly Raman, but valid for simulation)
                peakOmegaEff = m_config.m_peakRabiFrequency;
            }
            else
            {
                peakOmegaEff = (m_config.m_peakRabiFrequency * m_config.m_peakRabiFrequency)
                    / (2.0 * std::abs(m_config.m_commonDetuning));
            }

            // Total Area (Theta) = ∫ Ω_eff(t) dt = peakOmegaEff * σ * sqrt(π / 2)
            m_Sigma = m_config.m_targetTheta / (peakOmegaEff * std::sqrt(ConstexprMath::Pi / 2.0));

            m_tP = 4.0 * m_Sigma;
            m_tS = 4.0 * m_Sigma;

            // Full pulse duration (no truncation needed for Raman)
            m_TimeLimit = 8.0 * m_Sigma;
        }

        void setupStirapParameters()
        {
            // Adiabatic condition: Ω0 * σ >> 1. We use a heuristic for robust transfer.
            m_Sigma = 10.0 / m_config.m_peakRabiFrequency;

            if (m_config.m_protocol == TwoPhotonPulseProtocol::STIRAP)
            {
                m_tS = 3.0 * m_Sigma; // Stokes first
                m_tP = 5.0 * m_Sigma; // Pump second
            }
            else // Inverted STIRAP
            {
                m_tP = 3.0 * m_Sigma;
                m_tS = 5.0 * m_Sigma;
            }

            // If target is Pi, we run the full sequence.
            if (m_config.m_targetTheta >= ConstexprMath::Pi - 1e-7)
            {
                m_TimeLimit = std::max(m_tP, m_tS) + 4.0 * m_Sigma;
            }
            else
            {
                // Fractional STIRAP: stop when tan(theta/2) = Ωp(t)/Ωs(t)
                real_t ratio = std::tan(m_config.m_targetTheta / 2.0);
                real_t deltaT = m_tP - m_tS;

                // Analytical solution for t where Ωp(t)/Ωs(t) = ratio
                m_TimeLimit = ((m_tP + m_tS) / 2.0) + (m_Sigma * m_Sigma * std::log(ratio)) / (2.0 * deltaT);
            }
        }

    public:
        /**
            * @brief Returns the laser pulses at a given time.
            * @details Uses clamping at TimeLimit to maintain the mixing angle for fractional STIRAP.
            */
        std::tuple<LaserPulse, LaserPulse> operator()(real_t time) const noexcept
        {
            // Clamping ensures that if the solver exceeds TimeLimit, the state remains coherent.
            const real_t evalTime = std::min(time, m_TimeLimit);

            const real_t ampP = m_config.m_peakRabiFrequency * gaussian(evalTime, m_tP, m_Sigma);
            const real_t ampS = m_config.m_peakRabiFrequency * gaussian(evalTime, m_tS, m_Sigma);

            return {
                LaserPulse{ m_omegaP, ampP, m_config.m_pumpPhase },
                LaserPulse{ m_omegaS, ampS, 0.0 }
            };
        }

        /**
            * @brief The duration of the simulation required to achieve targetTheta.
            */
        constexpr real_t getTransitionTimeLimit() const noexcept
        {
            return m_TimeLimit;
        }
    };

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
        constexpr STIRAPLaser(const STIRAPConfig& config, real_t c1, real_t c2)
            : m_config(config)
        {
            // Target Rabi frequency parameters for defining the pulse area.
            m_Omega0 = Units::omegaAuFromHz(config.m_rabiFrequency_Hz);
            m_Sigma = 150.0 / m_Omega0;

            // Center-times for the delayed pulse sequence (Stokes precedes Pump).
            m_StokesTimeCenter = c1 * m_Sigma;
            m_PumpTimeCenter = c2 * m_Sigma;
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
