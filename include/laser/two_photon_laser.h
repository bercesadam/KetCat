#pragma once
#include <tuple>

#include "atomic_units.h"
#include "laser/laser_pulse.h"
#include "quantum_bit/neutral_atom_manifold.h"

namespace KetCat
{
    /// @brief Supported coherent two-photon control protocols.
    ///
    /// @details
    ///   These protocols generate time-dependent pump and Stokes laser
    ///   pulses for driving population transfer or coherent rotations
    ///   in a 3-level ladder system:
    ///
    ///     |0⟩ ↔ |1⟩ ↔ |2⟩
    ///
    ///   Protocol overview:
    ///
    ///     • STIRAP
    ///         Counter-intuitive pulse ordering:
    ///             Stokes → Pump
    ///
    ///         Adiabatic population transfer:
    ///             |0⟩ → |2⟩
    ///
    ///     • InvertedSTIRAP
    ///         Reverse pulse ordering:
    ///             Pump → Stokes
    ///
    ///         Adiabatic population transfer:
    ///             |2⟩ → |0⟩
    ///
    ///     • Simultaneous
    ///         Pump and Stokes overlap simultaneously,
    ///         producing an effective Raman coupling:
    ///
    ///             Ω_eff ≈ (Ωp Ωs) / (2Δ)
    enum class TwoPhotonProtocol
    {
        STIRAP,
        InvertedSTIRAP,
        Simultaneous
    };

    /// @brief Physical configuration for a coherent two-photon transition.
    ///
    /// @details
    ///   Defines the laser and atomic parameters required for generating
    ///   pump and Stokes pulse sequences.
    ///
    ///   Ladder structure:
    ///
    ///       |2⟩  ← Ωs
    ///        │
    ///       |1⟩
    ///        │
    ///       |0⟩  ← Ωp
    ///
    ///   Two-photon resonance condition:
    ///
    ///       (E₂ − E₀) = ℏ(ωp + ωs)
    ///
    ///   The common detuning Δ offsets the intermediate state:
    ///
    ///       Δ = (E₁ − E₀) − ℏωp
    struct TwoPhotonConfig
    {
        /// @brief Energy of logical ground state |0⟩.
        real_t m_Level1Energy;

        /// @brief Energy of intermediate state |1⟩.
        real_t m_Level2Energy;

        /// @brief Energy of target state |2⟩.
        real_t m_Level3Energy;

        /// @brief Transition dipole μ₀₁ = ⟨0|μ|1⟩.
        real_t m_Mu12;

        /// @brief Transition dipole μ₁₂ = ⟨1|μ|2⟩.
        real_t m_Mu23;

        /// @brief Peak single-photon Rabi frequency Ω₀.
        ///
        /// @details
        ///   Used as the maximum amplitude of both Gaussian pulses.
        real_t m_peakRabiFrequency;

        /// @brief Relative phase φ of the pump pulse.
        ///
        /// @details
        ///   Controls the effective Raman rotation axis.
        real_t m_pumpPhase;

        /// @brief Common detuning Δ from the intermediate state.
        ///
        /// @details
        ///   Positive detuning suppresses intermediate-state population
        ///   and produces an effective Raman coupling.
        real_t m_commonDetuning;

        /// @brief Target Bloch rotation angle θ.
        ///
        /// @details
        ///   Examples:
        ///
        ///     θ = π     → full population transfer
        ///     θ = π/2   → equal superposition
        real_t m_targetTheta;

        /// @brief Selected pulse protocol.
        TwoPhotonProtocol m_protocol = TwoPhotonProtocol::STIRAP;
    };

    /// @brief Time-dependent two-photon laser pulse generator.
    ///
    /// @details
    ///   Generates coherent pump and Stokes Gaussian laser pulses
    ///   for driving ladder-system dynamics.
    ///
    ///   Supported modes:
    ///
    ///     • STIRAP / Inverted STIRAP
    ///         Adiabatic population transfer
    ///
    ///     • Simultaneous Raman
    ///         Effective two-level coherent rotations
    ///
    ///   The generated pulses are intended for direct use with TDSE
    ///   propagators and RWA Hamiltonian builders.
    class TwoPhotonLaser
    {
        /// @brief Physical configuration of the pulse sequence.
        TwoPhotonConfig m_config;

        /// @brief Gaussian pulse width σ.
        real_t m_Sigma;

        /// @brief Temporal center of pump pulse.
        real_t m_tP;

        /// @brief Temporal center of Stokes pulse.
        real_t m_tS;

        /// @brief Total evolution duration.
        real_t m_TimeLimit;

		/// @brief Total time required for full population transfer.
        real_t m_FullTransferTime;

        /// @brief Pump laser angular frequency ωp.
        real_t m_omegaP;

        /// @brief Stokes laser angular frequency ωs.
        real_t m_omegaS;

        /// @brief Evaluate normalized Gaussian envelope.
        ///
        /// @param t
        ///   Evaluation time.
        ///
        /// @param t0
        ///   Pulse center.
        ///
        /// @param sigma
        ///   Gaussian width.
        ///
        /// @return
        ///   exp(-(t-t₀)² / σ²)
        static constexpr real_t gaussian(real_t t, real_t t0, real_t sigma) noexcept
        {
            if (sigma < 1e-12)
            {
                return 0.0;
            }

            const real_t x = (t - t0) / sigma;

            return std::exp(-x * x);
        }

    public:
        /// @brief Construct a two-photon pulse generator.
        ///
        /// @param config
        ///   Physical pulse configuration.
        ///
        /// @details
        ///   Laser frequencies are initialized to satisfy the
        ///   two-photon resonance condition:
        ///
        ///       (E₂ − E₀) = ℏ(ωp + ωs)
        ///
        ///   with opposite detuning applied to the two fields:
        ///
        ///       ωp = ω₀₁ + Δ
        ///       ωs = ω₁₂ − Δ
        constexpr TwoPhotonLaser(const TwoPhotonConfig& config)
            : m_config(config)
        {
            m_omegaP =
                (m_config.m_Level2Energy - m_config.m_Level1Energy)
                + m_config.m_commonDetuning;

            m_omegaS =
                (m_config.m_Level3Energy - m_config.m_Level2Energy)
                - m_config.m_commonDetuning;

            if (m_config.m_protocol == TwoPhotonProtocol::Simultaneous)
            {
                setupRamanParameters();
            }
            else
            {
                setupStirapParameters();
            }
        }

    private:
        /// @brief Configure simultaneous Raman pulse parameters.
        ///
        /// @details
        ///   The effective Raman coupling is approximated as:
        ///
        ///       Ω_eff(t) = Ωp(t) Ωs(t) / (2Δ)
        ///
        ///   For identical overlapping Gaussian pulses:
        ///
        ///       Ω_eff(t) ∝ exp(-2(t-t₀)² / σ²)
        ///
        ///   The pulse width σ is selected such that:
        ///
        ///       θ = ∫ Ω_eff(t) dt
        ///
        ///   where θ is the target Bloch rotation angle.
        void setupRamanParameters()
        {
            real_t peakOmegaEff;

            if (std::abs(m_config.m_commonDetuning) < 1e-10)
            {
                /// Resonant limit
                peakOmegaEff = m_config.m_peakRabiFrequency;
            }
            else
            {
                peakOmegaEff =
                    (m_config.m_peakRabiFrequency
                        * m_config.m_peakRabiFrequency)
                    / (2.0 * std::abs(m_config.m_commonDetuning));
            }

            /// θ = Ω_eff σ √(π/2)
            m_Sigma =
                m_config.m_targetTheta
                / (peakOmegaEff * std::sqrt(ConstexprMath::Pi / 2.0));

            m_tP = 4.0 * m_Sigma;
            m_tS = 4.0 * m_Sigma;

            m_TimeLimit = 8.0 * m_Sigma;
        }

        /// @brief Configure STIRAP pulse timing parameters.
        ///
        /// @details
        ///   The pulse width σ is selected heuristically to satisfy
        ///   the adiabatic condition:
        ///
        ///       Ω₀ σ ≫ 1
        ///
        ///   Pulse ordering:
        ///
        ///     • STIRAP:
        ///         Stokes before Pump
        ///
        ///     • Inverted STIRAP:
        ///         Pump before Stokes
        ///
        ///   For fractional transfer:
        ///
        ///       tan(θ/2) = Ωp(t) / Ωs(t)
        ///
        ///   the sequence terminates early at the corresponding
        ///   mixing angle.
        void setupStirapParameters()
        {
            m_Sigma = 10.0 / m_config.m_peakRabiFrequency;

            if (m_config.m_protocol == TwoPhotonProtocol::STIRAP)
            {
                m_tS = 2.5 * m_Sigma;
                m_tP = 4.5 * m_Sigma;
            }
            else
            {
                m_tS = 4.5 * m_Sigma;
                m_tP = 2.5 * m_Sigma;
            }

            m_FullTransferTime =
                std::max(m_tP, m_tS)
                + 3.0 * m_Sigma;

            if (ConstexprMath::floatNear(m_config.m_targetTheta, ConstexprMath::Pi))
            {
				m_TimeLimit = m_FullTransferTime;
            }
            else
            {
                const real_t ratio =
                    std::tan(m_config.m_targetTheta / 2.0);

                const real_t deltaT = m_tP - m_tS;

                m_TimeLimit =
                    ((m_tP + m_tS) / 2.0)
                    + (m_Sigma * m_Sigma * std::log(ratio))
                    / (2.0 * deltaT);
            }
        }

    public:
        /// @brief Evaluate pump and Stokes pulses at time t.
        ///
        /// @param time
        ///   Simulation time.
        ///
        /// @return
        ///   Tuple:
        ///
        ///     (PumpLaser, StokesLaser)
        ///
        /// @details
        ///   Time evaluation is clamped at m_TimeLimit.
        ///
        ///   This preserves the final mixing angle for fractional
        ///   STIRAP operations if the solver continues evolving
        ///   beyond the designed pulse duration.
        std::tuple<LaserPulse, LaserPulse> operator()(real_t time) const noexcept
        {
            const real_t evalTime =
                std::min(time, m_TimeLimit);

            const real_t ampP =
                m_config.m_peakRabiFrequency
                * gaussian(evalTime, m_tP, m_Sigma);

            const real_t ampS =
                m_config.m_peakRabiFrequency
                * gaussian(evalTime, m_tS, m_Sigma);

            return {
                LaserPulse{ m_omegaP, ampP, m_config.m_pumpPhase },
                LaserPulse{ m_omegaS, ampS, 0.0 }
            };
        }

        /// @brief Get total pulse duration required for target transfer.
        ///
        /// @return
        ///   Final transition time limit.
        constexpr real_t getTransitionTimeLimit() const noexcept
        {
            return m_TimeLimit;
        }

		/// @brief Get total time required for full population transfer.
        ///
		/// @return 
		///   Time required to achieve complete population transfer (θ = π) for STIRAP protocols.
        constexpr real_t getFullTransferTime() const noexcept
        {
            return m_FullTransferTime;
		}
    };
}