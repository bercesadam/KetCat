#pragma once
#include <tuple>
#include <cmath>
#include <algorithm>

#include "atomic_units.h"
#include "laser/laser_pulse.h"
#include "operation_space/neutral_atom_manifold.h"


namespace KetCat
{
    /// @brief Supported coherent two-photon control protocols.
    enum class TwoPhotonProtocol
    {
        STIRAP,
        InvertedSTIRAP,
        Simultaneous
    };

    /// @brief Physical configuration for a coherent two-photon transition.
    struct TwoPhotonConfig
    {
        real_t m_Level1Energy;        ///< Energy of state |1⟩ (ground state, logical |0⟩)
        real_t m_Level2Energy;        ///< Energy of state |2⟩ (intermediate state)
        real_t m_Level3Energy;        ///< Energy of state |3⟩ (excited state, logical |1⟩)
        real_t m_LevelRydbergEnergy;  ///< Energy of the Rydberg state
        real_t m_Mu12;                ///< Dipole μ₁₂
        real_t m_Mu23;                ///< Dipole μ₂₃
        real_t m_Mu3r;                ///< Dipole μ₃ᵣ
        real_t m_peakRabiFrequency;   ///< Max Ω₀
        real_t m_pumpPhase;           ///< Phase φ
        real_t m_commonDetuning;      ///< Detuning Δ
        real_t m_targetTheta;         ///< Rotation angle θ
        TwoPhotonProtocol m_protocol = TwoPhotonProtocol::STIRAP;
    };

    /// @brief Pre-calculated parameters for pulse evaluation.
    /// @details Contains the precalculated values required for the Gaussian envelope generation.
    struct EnvelopeParams
    {
        real_t m_omegaP;    ///< Pump frequency
        real_t m_omegaS;    ///< Stokes frequency
        real_t m_peakRabi;  ///< Peak amplitude
        real_t m_phaseP;    ///< Pump phase
        real_t m_sigma;     ///< Gaussian width
        real_t m_tP;        ///< Pump center
        real_t m_tS;        ///< Stokes center
        real_t m_tLimit;    ///< Evaluation cutoff
        real_t m_tStart;    ///< Start time in case of fractional STIRAP
    };

    /// @brief Lightweight time-dependent two-photon laser pulse evaluator.
    ///
    /// @details
    ///    This class is intended to be distributed to solvers. It contains 
    ///    no configuration logic, only the mathematical evaluation of the 
    ///    pre-calculated pulse sequence.
    class TwoPhotonLaserEnvelope
    {
        EnvelopeParams m_Parameters;

        /// @brief Evaluate normalized Gaussian envelope.
        static constexpr real_t gaussian(real_t t, real_t t0, real_t sigma) noexcept
        {
            if (sigma < 1e-12)
            {
                return 0.0;
            }
            const real_t x = (t - t0) / sigma;
            return ConstexprMath::exp(-x * x);
        }

    public:
        /// @brief Construct from pre-calculated parameters.
        constexpr TwoPhotonLaserEnvelope(const EnvelopeParams& params) noexcept
            : m_Parameters(params) {
        }

        /// @brief Evaluate pump and Stokes pulses at time t.
        /// @return Tuple: (PumpLaser, StokesLaser)
        std::tuple<LaserPulse, LaserPulse> operator()(real_t time) const noexcept
        {
            const real_t evalTime = std::min(time, m_Parameters.m_tLimit);

            const real_t ampP = m_Parameters.m_peakRabi *
                gaussian(evalTime, m_Parameters.m_tP, m_Parameters.m_sigma);
            const real_t ampS = m_Parameters.m_peakRabi *
                gaussian(evalTime, m_Parameters.m_tS, m_Parameters.m_sigma);

            return {
                LaserPulse{ m_Parameters.m_omegaP, ampP, m_Parameters.m_phaseP },
                LaserPulse{ m_Parameters.m_omegaS, ampS, 0.0 }
            };
        }

        /// @brief Get total pulse duration required for target transfer.
        constexpr real_t getTransitionTimeLimit() const noexcept { return m_Parameters.m_tLimit; }

        real_t setStartTime(const real_t startTime) noexcept { m_Parameters.m_tStart = startTime; }
        constexpr real_t getStartTime() const noexcept { return m_Parameters.m_tStart; }
    };

    /// @brief Builder class for TwoPhotonLaserEnvelope.
    ///
    /// @details
    ///    Encapsulates the heuristic logic for setting up pulse widths,
    ///    timings, and adiabatic conditions based on physical configuration.
    class TwoPhotonPulseBuilder
    {
        TwoPhotonConfig m_config;

    public:
        explicit TwoPhotonPulseBuilder(const TwoPhotonConfig& config)
            : m_config(config) {
        }

        /// @brief Compute parameters and generate the final envelope evaluator.
        TwoPhotonLaserEnvelope build() const
        {
            EnvelopeParams Parameters;

            // 1. Calculate base frequencies
            Parameters.m_omegaP = (m_config.m_Level2Energy - m_config.m_Level1Energy) + m_config.m_commonDetuning;
            Parameters.m_omegaS = (m_config.m_Level3Energy - m_config.m_Level2Energy) - m_config.m_commonDetuning;

            Parameters.m_peakRabi = m_config.m_peakRabiFrequency;
            Parameters.m_phaseP = m_config.m_pumpPhase;

            // 2. Protocol-specific timing calculation
            if (m_config.m_protocol == TwoPhotonProtocol::Simultaneous)
            {
                setupRaman(Parameters);
            }
            else
            {
                setupStirap(Parameters);
            }

            return TwoPhotonLaserEnvelope(Parameters);
        }

    private:
        // Curently not used but kept for possible future uses
        void setupRaman(EnvelopeParams& Parameters) const
        {
            real_t peakOmegaEff;
            if (ConstexprMath::abs(m_config.m_commonDetuning) < 1e-10)
            {
                peakOmegaEff = m_config.m_peakRabiFrequency;
            }
            else
            {
                peakOmegaEff = (m_config.m_peakRabiFrequency * m_config.m_peakRabiFrequency)
                    / (2.0 * ConstexprMath::abs(m_config.m_commonDetuning));
            }

            Parameters.m_sigma = m_config.m_targetTheta / (peakOmegaEff * ConstexprMath::sqrt(ConstexprMath::Pi / 2.0));
            Parameters.m_tP = 4.0 * Parameters.m_sigma;
            Parameters.m_tS = 4.0 * Parameters.m_sigma;
            Parameters.m_tLimit = 8.0 * Parameters.m_sigma;
        }

        void setupStirap(EnvelopeParams& Parameters) const
        {
            Parameters.m_sigma = 10.0 / m_config.m_peakRabiFrequency;

            if (m_config.m_protocol == TwoPhotonProtocol::STIRAP)
            {
                Parameters.m_tS = 2.5 * Parameters.m_sigma;
                Parameters.m_tP = 4.5 * Parameters.m_sigma;
            }
            else
            {
                Parameters.m_tS = 4.5 * Parameters.m_sigma;
                Parameters.m_tP = 2.5 * Parameters.m_sigma;
            }

            const real_t fullTransferTime =
                std::max(Parameters.m_tP, Parameters.m_tS) + 3.0 * Parameters.m_sigma;

            if (ConstexprMath::floatNear(m_config.m_targetTheta, ConstexprMath::Pi))
            {
                Parameters.m_tLimit = fullTransferTime;
            }
            else
            {
                // This is an empirical safety margin which lets the laser run for 
                // a little bit longer than the time yielded from the theoretical formula
                // This is a simple solution to avoid issues because of abrupted pulses
                // a better solution would be a gentle logarithmic tail to be developed later
                constexpr real_t SafetyMarginRatio = 0.005;

                const real_t Ratio = ConstexprMath::tan(m_config.m_targetTheta / 2.0);
                const real_t DeltaT = Parameters.m_tP - Parameters.m_tS;

                Parameters.m_tLimit = ((Parameters.m_tP + Parameters.m_tS) / 2.0)
                    + (Parameters.m_sigma * Parameters.m_sigma * ConstexprMath::log(Ratio)) / (2.0 * DeltaT)
                    + (fullTransferTime * SafetyMarginRatio);
            }
        }
    };
}