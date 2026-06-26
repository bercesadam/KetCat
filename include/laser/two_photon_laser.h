#pragma once
#include <tuple>
#include <cmath>
#include <algorithm>

#include "laser/laser_pulse.h"


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
        natural_t m_GroundLevelIndex; ///< The first affected bases state index to fill the resulting laser array
        real_t m_Level1Energy;        ///< Energy of state |1⟩ (ground state, logical |0⟩)
        real_t m_Level2Energy;        ///< Energy of state |2⟩ (intermediate state)
        real_t m_Level3Energy;        ///< Energy of state |3⟩ (excited state, logical |1⟩)
        real_t m_LevelRydbergEnergy;  ///< Energy of the Rydberg state
        real_t m_Mu12;                ///< Dipole μ₁₂
        real_t m_Mu23;                ///< Dipole μ₂₃
        real_t m_Mu3r;                ///< Dipole μ₃ᵣ
        real_t m_peakRabiFrequencyP;   ///< Max Ω₀ (Pump)
        real_t m_peakRabiFrequencyS;   ///< Max Ω₀ (Stokes)
        real_t m_pumpPhase;           ///< Phase φ
        real_t m_pumpPhase2;
        real_t m_pumpPhase2Timing;
        real_t m_commonDetuning;      ///< Detuning Δ
        real_t m_targetTheta;         ///< Rotation angle θ
        TwoPhotonProtocol m_protocol = TwoPhotonProtocol::STIRAP;
    };

    /// @brief Pre-calculated parameters for pulse evaluation.
    /// @details Contains the precalculated values required for the Gaussian envelope generation.
    namespace detail
    {
        struct EnvelopeParams
        {
            natural_t m_GroundLevelIndex;  ///< The first affected bases state index to fill the resulting laser array
            real_t m_omegaP;         ///< Pump frequency
            real_t m_omegaS;         ///< Stokes frequency
            real_t m_peakRabiP;       ///< Peak amplitude for the pump laser
            real_t m_peakRabiS;       ///< Peak amplitude for the Stokes laser
            real_t m_phaseP;         ///< Pump phase
            real_t m_phaseP2;
            real_t m_phaseP2Time;
            real_t m_sigma;          ///< Gaussian width
            real_t m_tP;             ///< Pump center
            real_t m_tS;             ///< Stokes center
            real_t m_tLimit;         ///< Evaluation cutoff
            real_t m_tStart;         ///< Start time in case of fractional STIRAP
            real_t m_PiTransferTime; ///< Theoretical transfer time needed for a full Pi transition
        };
    }

    /// @brief Lightweight time-dependent two-photon laser pulse evaluator.
    ///
    /// @details
    ///    This class is intended to be distributed to solvers. It contains 
    ///    no configuration logic, only the mathematical evaluation of the 
    ///    pre-calculated pulse sequence.
    class TwoPhotonLaserEnvelope
    {
        detail::EnvelopeParams m_Parameters;

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
        constexpr TwoPhotonLaserEnvelope(const detail::EnvelopeParams& params) noexcept
            : m_Parameters(params) {
        }

        constexpr natural_t getGroundLevelIndex() const noexcept
        {
            return m_Parameters.m_GroundLevelIndex;
        }

        /// @brief Evaluate pump and Stokes pulses at time t.
        /// @return Tuple: (PumpLaser, StokesLaser)
        TwoPhotonDrive operator()(real_t time) const noexcept
        {
            const real_t AmpP = m_Parameters.m_peakRabiP *
                gaussian(time, m_Parameters.m_tP, m_Parameters.m_sigma);
            const real_t AmpS = m_Parameters.m_peakRabiS *
                gaussian(time, m_Parameters.m_tS, m_Parameters.m_sigma);

            real_t ActualPhase = m_Parameters.m_phaseP;
            if ((time / m_Parameters.m_tLimit) >= m_Parameters.m_phaseP2Time)
            {
                ActualPhase = m_Parameters.m_phaseP2;
            }
           

            return
            {
                m_Parameters.m_GroundLevelIndex,
                LaserPulse{ m_Parameters.m_omegaP, AmpP, ActualPhase },
                LaserPulse{ m_Parameters.m_omegaS, AmpS, 0.0 }
            };
        }

        /// @brief Get total pulse duration required for target transfer.
        constexpr real_t getTransitionTimeLimit() const noexcept { return m_Parameters.m_tLimit; }

        /// @brief Get total pulse duration required for a theoretical (or actual) full Pi transfer.
        constexpr real_t getFullTransferTime() const noexcept { return m_Parameters.m_PiTransferTime; }

        constexpr void setStartTime(const real_t startTime) noexcept { m_Parameters.m_tStart = startTime; }
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
            detail::EnvelopeParams Parameters;

			// Copy the index of the first affected base state which will be used to fill the lasers array - TODO move to laser_pulse type
			Parameters.m_GroundLevelIndex = m_config.m_GroundLevelIndex;

            // 1. Calculate base frequencies
            Parameters.m_omegaP = (m_config.m_Level2Energy - m_config.m_Level1Energy) + m_config.m_commonDetuning;
            Parameters.m_omegaS = (m_config.m_Level3Energy - m_config.m_Level2Energy) - m_config.m_commonDetuning;

			// 2. Calculate peak Rabi frequencies based on the provided peak Rabi frequency and the dipole moments.
            Parameters.m_peakRabiP = m_config.m_peakRabiFrequencyP;
            Parameters.m_peakRabiS = m_config.m_peakRabiFrequencyS;

			// 3. Copy the pump phase from the configuration, which already includes the rotating frame correction.
            Parameters.m_phaseP = m_config.m_pumpPhase;
            Parameters.m_phaseP2 = m_config.m_pumpPhase2;
            Parameters.m_phaseP2Time = m_config.m_pumpPhase2Timing;

            // 4. Protocol-specific timing calculation
            if (m_config.m_protocol == TwoPhotonProtocol::Simultaneous)
            {
                setupRaman(Parameters);
            }
            else
            {
            }

            return TwoPhotonLaserEnvelope(Parameters);
        }

    private:
        void setupRaman(detail::EnvelopeParams& Parameters) const
        {
            real_t geomPeakRabi = std::sqrt(Parameters.m_peakRabiP * Parameters.m_peakRabiS);

            if (std::abs(m_config.m_commonDetuning) < 1e-10)
            {
				// Resonant Raman case: the effective Rabi frequency is approximately the geometric mean of the two peaks.
                real_t peakOmegaEff = geomPeakRabi * ConstexprMath::sqrt(2.0);
                Parameters.m_sigma = m_config.m_targetTheta / (peakOmegaEff * ConstexprMath::sqrt(ConstexprMath::Pi / 2.0));
            }
            else
            {
				// Off-resonant Raman case: the effective Rabi frequency is reduced by the detuning, leading to a longer required pulse duration.
                real_t peakOmegaEff = (geomPeakRabi * geomPeakRabi) / (2.0 * ConstexprMath::abs(m_config.m_commonDetuning));

				// The pulse width (sigma) is inversely proportional to the effective Rabi frequency, which is reduced in the off-resonant case.
                Parameters.m_sigma = m_config.m_targetTheta / (peakOmegaEff * ConstexprMath::sqrt(ConstexprMath::Pi / 2.0));
                // To ensure sufficient adiabaticity in the off-resonant case, we apply an additional safety factor to the pulse width.
                Parameters.m_sigma *= ConstexprMath::sqrt(2.0);
            }

            Parameters.m_tP = 4.0 * Parameters.m_sigma;
            Parameters.m_tS = 4.0 * Parameters.m_sigma;
            Parameters.m_tLimit = 8.0 * Parameters.m_sigma;
        }
    };
}