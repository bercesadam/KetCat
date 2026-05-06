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
}
