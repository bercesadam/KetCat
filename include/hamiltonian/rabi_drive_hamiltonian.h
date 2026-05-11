#pragma once 
#include "core_types.h"
#include "laser/laser_pulse.h"


namespace KetCat
{
    /// @brief Multi-laser tridiagonal Hamiltonian generator in the Rotating Wave Approximation (RWA).
    ///
    /// @details
    ///   This class constructs the effective Hamiltonian for an N-level ladder system
    ///   driven by multiple coherent laser fields (e.g. |0⟩ → |1⟩ → |2⟩).
    ///
    ///   The system is transformed into a rotating frame, where fast oscillating
    ///   terms are removed and the dynamics becomes slowly varying.
    ///
    ///   Matrix structure (tridiagonal):
    ///
    ///     • Main diagonal:
    ///         Hᵢᵢ = Δᵢ + Σ ΔE_stark
    ///
    ///     • Off-diagonals:
    ///         Hᵢ,ᵢ₊₁ = Ωᵢ / 2
    ///         Hᵢ₊₁,ᵢ = (Ωᵢ / 2)*
    ///
    ///   where:
    ///
    ///     Δᵢ        = cumulative multi-photon detuning
    ///     Ωᵢ        = Rabi frequency of transition |i⟩ ↔ |i+1⟩
    ///     ΔE_stark  = AC Stark shift from off-resonant couplings
    ///
    ///   Key properties:
    ///
    ///     • Only nearest-neighbor transitions are explicitly coupled
    ///     • Arbitrary off-resonant couplings contribute via Stark shifts
    ///     • Hamiltonian remains Hermitian by construction
    ///
    /// @tparam LevelCount
    ///   Total number of energy levels in the simulation manifold.
    template<natural_t LevelCount>
    class MultiRwaRabiHamiltonian
    {
    public:
        using ReducedMatrix = tridiagonal_matrix_t<LevelCount>;
        using FullDipoleMatrix = matrix_t<LevelCount>;

    private:
        /// @brief Bare energy levels of the system (e.g. Hartree units).
        std::array<real_t, LevelCount> m_Energies{};

        /// @brief Full dipole transition matrix μᵢⱼ.
        FullDipoleMatrix m_DipoleMatrix{};

        /// @brief Internal tridiagonal Hamiltonian representation.
        tridiagonal_matrix_t<LevelCount> m_hamiltonianMatrix;

    public:
        /// @brief Construct and initialize the RWA Hamiltonian.
        ///
        /// @param energies
        ///   Bare energy levels of the system.
        ///
        /// @param dipoleMatrix
        ///   Transition dipole matrix μᵢⱼ.
        ///
        /// @param lasers
        ///   Array of laser pulses, where lasers[i] drives |i⟩ → |i+1⟩.
        ///
        /// @details
        ///   Initializes both diagonal (detuning + Stark shifts) and
        ///   off-diagonal (Rabi couplings) terms.
        constexpr MultiRwaRabiHamiltonian(
            const std::array<real_t, LevelCount>& energies,
            const FullDipoleMatrix& dipoleMatrix,
            const std::array<LaserPulse, LevelCount - 1>& lasers) noexcept
            : m_Energies(energies),
            m_DipoleMatrix(dipoleMatrix)
        {
            m_hamiltonianMatrix = {};
            updateMainDiagonal(lasers);
            updateOffDiagonal(lasers);
        }

        /// @brief Get the current Hamiltonian matrix.
        ///
        /// @return
        ///   Reference to the internal tridiagonal matrix.
        constexpr const tridiagonal_matrix_t<LevelCount>& getMatrix() const noexcept
        {
            return m_hamiltonianMatrix;
        }

        //private:
            /// @brief Update diagonal elements: detuning + AC Stark shifts.
            ///
            /// @param lasers
            ///   Current laser configuration.
            ///
            /// @details
            ///   The diagonal term for level i is:
            ///
            ///     Hᵢᵢ = Δᵢ + ΔE_stark
            ///
            ///   where:
            ///
            ///     Δᵢ = (Eᵢ − E₀) − Σ ωⱼ
            ///
            ///   is the cumulative detuning in the rotating frame.
            ///
            ///   The AC Stark shift is computed using second-order perturbation:
            ///
            ///     ΔE_stark ≈ Σ |Ωᵢⱼ/2|² / Δ
            ///
            ///   including:
            ///     • Rotating (resonant) contributions
            ///     • Counter-rotating (Bloch–Siegert) terms
        constexpr void updateMainDiagonal(const std::array<LaserPulse, LevelCount - 1>& lasers) noexcept
        {
            real_t CumulativeOmega = 0.0;

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                /// Accumulate rotating-frame frequency shift
                if (i > 0)
                {
                    CumulativeOmega += lasers[i - 1].m_omega;
                }

                const real_t RelativeEnergy = m_Energies[i] - m_Energies[0];
                const real_t Detuning = RelativeEnergy - CumulativeOmega;

                /// --- AC Stark shift calculation ---
                real_t StarkShift = 0.0;

                for (natural_t level = 0; level < LevelCount - 1; ++level)
                {
                    const auto& Laser = lasers[level];

                    for (natural_t j = 0; j < LevelCount; ++j)
                    {
                        if (i == j)
                        {
                            continue;
                        }

                        const real_t DipoleElm = m_DipoleMatrix[i][j].re;

                        /// Half Rabi frequency:
                        ///   Ωᵢⱼ/2 = −½ · μᵢⱼ · E
                        const real_t OmegaRabiHalf = -0.5 * DipoleElm * Laser.m_amplitude;

                        const real_t DeltaE = m_Energies[i] - m_Energies[j];

                        /// Rotating (near-resonant) contribution
                        const real_t ResonantDenom = DeltaE - Laser.m_omega;
                        if (ConstexprMath::abs(ResonantDenom) > 1e-9)
                        {
                            StarkShift += (OmegaRabiHalf * OmegaRabiHalf) / ResonantDenom;
                        }

                        /// Counter-rotating (Bloch–Siegert) contribution
                        const real_t anti_ResonantDenom = DeltaE + Laser.m_omega;
                        if (ConstexprMath::abs(anti_ResonantDenom) > 1e-9)
                        {
                            StarkShift += (OmegaRabiHalf * OmegaRabiHalf) / anti_ResonantDenom;
                        }
                    }
                }

                m_hamiltonianMatrix[MAINDIAGONAL][i] =
                    complex_t::fromReal(Detuning + StarkShift);
            }
        }

        /// @brief Update off-diagonal elements (coherent couplings).
        ///
        /// @param lasers
        ///   Current laser configuration.
        ///
        /// @details
        ///   Each laser drives a transition:
        ///
        ///     Hᵢ,ᵢ₊₁ = −½ · μᵢ,ᵢ₊₁ · E · e^{iφ}
        ///
        ///   The Hermitian conjugate ensures:
        ///
        ///     Hᵢ₊₁,ᵢ = (Hᵢ,ᵢ₊₁)*
        void updateOffDiagonal(const std::array<LaserPulse, LevelCount - 1>& lasers)
        {
            for (natural_t i = 0; i < LevelCount - 1; ++i)
            {
                const LaserPulse& Laser = lasers[i];

                /// Complex phase factor e^{iφ}
                const complex_t PhaseFactor =
                {
                    std::cos(Laser.m_phase),
                    -std::sin(Laser.m_phase)
                };

                /// Coupling term:
                ///   Ω/2 = −½ · μ · E · e^{iφ}
                const complex_t Coupling = m_DipoleMatrix[i][i + 1]
                    * (0.5 * Laser.m_amplitude)
                    * PhaseFactor;

                m_hamiltonianMatrix[SUPERDIAGONAL][i] = Coupling;
                m_hamiltonianMatrix[SUBDIAGONAL][i + 1] = Coupling.conj();
            }
        }
    };
}