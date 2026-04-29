#pragma once 
#include "core_types.h"
#include "laser/laser.h"


namespace KetCat
{
    /// @brief Multi-laser Tridiagonal Hamiltonian generator in the Rotating Wave Approximation (RWA).
    ///
    /// @details
    /// This class constructs the Hamiltonian for an N-level "ladder" system (e.g., Cs 6s → 7p → 7s).
    /// The Hamiltonian is represented in a rotating frame, transforming the time-dependent
    /// Schrödinger equation into a time-independent (or slowly varying) eigenvalue problem.
    ///
    /// The matrix structure is tridiagonal:
    /// - **Main Diagonal:** Cumulative Detunings and second-order AC Stark shifts.
    /// - **Off-Diagonals:** Coherent Rabi couplings between adjacent levels.
    ///
    /// @tparam LevelCount Total number of energy levels in the simulation manifold.
    template<natural_t LevelCount>
    class MultiRwaRabiHamiltonian
    {
    public:
        using ReducedMatrix = tridiagonal_matrix_t<LevelCount>;
        using FullDipoleMatrix = matrix_t<LevelCount>;

    private:
        std::array<real_t, LevelCount> m_Energies{};
        FullDipoleMatrix m_DipoleMatrix{};
        std::array<LaserPulse, LevelCount - 1> m_Lasers{};

        tridiagonal_matrix_t<LevelCount> m_hamiltonianMatrix;

    public:
        /// @brief Constructs and initializes the RWA Hamiltonian.
        ///
        /// @param energies     Baseline energies of the states (e.g., Hartree).
        /// @param dipoleMatrix Matrix containing transition dipole moments μ_ij.
        /// @param lasers       Array of lasers where lasers[i] drives the |i⟩ → |i+1⟩ transition.
        constexpr MultiRwaRabiHamiltonian(
            const std::array<real_t, LevelCount>& energies,
            const FullDipoleMatrix& dipoleMatrix,
            const std::array<LaserPulse, LevelCount - 1>& lasers) noexcept
            : m_Energies(energies),
            m_DipoleMatrix(dipoleMatrix),
            m_Lasers(lasers)
        {
            calculateMatrix();
        }

        /// @return The constructed tridiagonal Hamiltonian matrix.
        constexpr tridiagonal_matrix_t<LevelCount> getMatrix() const noexcept
        {
            return m_hamiltonianMatrix;
        }

    private:
        /// @brief Computes the matrix elements, accounting for Detuning and AC Stark shifts.
        ///
        /// @details
        /// The diagonal elements are calculated as:
        ///     H_ii = Δ_i + Σ ΔE_Stark
        /// where Δ_i is the cumulative multi-photon Detuning.
        ///
        /// The off-diagonal elements (Rabi coupling) are:
        ///     Ω_i / 2 = -0.5 * μ_i,i+1 * ε * exp(iφ)
        constexpr void calculateMatrix() noexcept
        {
            m_hamiltonianMatrix = {};

            // --- 1. Diagonal Elements: Detuning + AC Stark-shift ---
            real_t CumulativeOmega = 0.0;

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                // To transform to the RWA frame, each level 'i' is shifted by the 
                // sum of the frequencies of the lasers used to reach it.
                if (i > 0)
                {
                    CumulativeOmega += m_Lasers[i - 1].m_omega;
                }

                const real_t RelativeEnergy = m_Energies[i] - m_Energies[0];
                const real_t Detuning = RelativeEnergy - CumulativeOmega;

                // --- AC Stark-shift calculation (Second-order Perturbation) ---
                // We account for non-resonant couplings from all lasers to all levels,
                // including the counter-rotating terms.
                real_t StarkShift = 0.0;
                for (natural_t level = 0; level < LevelCount - 1; ++level)
                {
                    const auto& Laser = m_Lasers[level];

                    for (natural_t j = 0; j < LevelCount; ++j)
                    {
						// Skip the diagonal term (i == j) since it does not contribute to the Stark shift.
                        if (i == j)
                        {
                            continue;
                        }

                        const real_t DipoleElm = m_DipoleMatrix[i][j].re;
                        const real_t OmegaRabiHalf = -0.5 * DipoleElm * Laser.m_amplitude;
                        const real_t DeltaE = m_Energies[i] - m_Energies[j];

                        // A. Resonant (Rotating) term:
                        // If the laser is resonant with the i-j transition, this term 
                        // is already handled by the off-diagonal coupling. We only 
                        // add the shift for non-resonant (far-detuned) transitions.
                        const real_t ResonantDenom = DeltaE - Laser.m_omega;
                        if (ConstexprMath::abs(ResonantDenom) > 1e-9)
                        {
                            StarkShift += (OmegaRabiHalf * OmegaRabiHalf) / ResonantDenom;
                        }

                        // B. Anti-resonant (Counter-rotating) term:
                        // This accounts for the Bloch-Siegert shift contribution.
                        const real_t anti_ResonantDenom = DeltaE + Laser.m_omega;
                        if (ConstexprMath::abs(anti_ResonantDenom) > 1e-9)
                        {
                            StarkShift += (OmegaRabiHalf * OmegaRabiHalf) / anti_ResonantDenom;
                        }
                    }
                }

                m_hamiltonianMatrix[MAINDIAGONAL][i] = complex_t::fromReal(Detuning + StarkShift);
            }

            // --- 2. Off-Diagonal Elements: Direct Rabi Couplings ---
            if constexpr (LevelCount > 1)
            {
                for (natural_t i = 0; i < LevelCount - 1; ++i)
                {
                    const LaserPulse& Laser = m_Lasers[i];

                    // Complex phase factor for the i -> i+1 driving field.
                    const complex_t PhaseFactor =
                    {
                        std::cos(Laser.m_phase),
                        std::sin(Laser.m_phase)
                    };

                    // The coupling term (half-Rabi frequency) including the transition phase.
                    // H_i,i+1 = -0.5 * d * E * e^(iφ)
                    const complex_t Coupling = m_DipoleMatrix[i][i + 1]
                        * (-0.5 * Laser.m_amplitude)
                        * PhaseFactor;

                    m_hamiltonianMatrix[SUPERDIAGONAL][i] = Coupling;
                    m_hamiltonianMatrix[SUBDIAGONAL][i + 1] = Coupling.conj(); // Hermiticity
                }
            }
        }
    };
}