#pragma once 
#include "core_types.h"

namespace KetCat
{
    /// @brief Tridiagonal Hamiltonian generator for Rabi drive in RWA frame
    ///
    /// @details
    ///   This class constructs a tridiagonal Hamiltonian matrix in the rotating frame.
    ///   It extracts only the nearest-neighbor couplings from a full dipole matrix to 
    ///   remain compatible with Thomas-algorithm (TDMA) based solvers.
    ///
    ///   H_tridiag = 
    ///   [ Δ₀   Ω₀₁  0   ]
    ///   [ Ω₁₀  Δ₁   Ω₁₂ ]
    ///   [ 0    Ω₂₁  Δ₂  ]
    ///
    /// @tparam LevelCount The number of levels in the reduced energy space.
    template<natural_t LevelCount>
    class RwaRabiHamiltonian
    {
    public:
        /// The reduced matrix type required by the Crank-Nicolson solver.
        using ReducedMatrix = tridiagonal_matrix_t<LevelCount>;
        using FullDipoleMatrix = matrix_t<LevelCount>;

    private:
        std::array<real_t, LevelCount> m_Energies{};
        FullDipoleMatrix m_DipoleMatrix{};
        real_t m_DriveOmega{};
        real_t m_DriveAmplitude{};
        real_t m_ReferenceLevel{};

        tridiagonal_matrix_t<LevelCount> m_hamiltonianMatrix;

    public:
        /// @brief Constructs the tridiagonal RWA generator.
        ///
        /// @param energies       Energies from ReducedEnergySpace::getEnergies().
        /// @param dipoleMatrix   Full dipole matrix from buildRadialDipoleMatrix().
        /// @param driveOmega     Laser angular frequency.
        /// @param referenceLevel The state index that defines the rotating frame.
        constexpr RwaRabiHamiltonian(
            const std::array<real_t, LevelCount>& energies,
            const FullDipoleMatrix& dipoleMatrix,
            const real_t driveOmega,
            const real_t driveAmplitude,
            const natural_t referenceLevel = 0) noexcept
            : m_Energies(energies),
            m_DipoleMatrix(dipoleMatrix),
            m_DriveOmega(driveOmega),
            m_DriveAmplitude(driveAmplitude),
            m_ReferenceLevel(referenceLevel)
        {
            calculateMatrix();
        }

        constexpr tridiagonal_matrix_t<LevelCount> getMatrix() const noexcept
        {
            return m_hamiltonianMatrix;
        }

    private:

        /// @brief Builds a tridiagonal Hamiltonian for the current time step.
        ///
        /// @tparam DriveFunctor  Functor for the electric field envelope ε(t).
        /// @param time           Current simulation time.
        /// @param drive          Field envelope function.
        constexpr void calculateMatrix() noexcept
        {
            m_hamiltonianMatrix = {};

            /// 1. Populate the Main Diagonal (Detunings + optional DC Stark shifts)
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                /// Energy in the rotating frame: Δᵢ = Eᵢ - (i - i_ref)·ω
                const real_t detuning =
                    m_Energies[i] - static_cast<real_t>(i - m_ReferenceLevel) * m_DriveOmega;

                /// If the diagonal dipole element μᵢᵢ is non-zero, include the DC Stark shift
                const real_t dcShift = 0.0; // -m_DipoleMatrix[i][i].re * m_DriveAmplitude;

                m_hamiltonianMatrix[MAINDIAGONAL][i] = complex_t::fromReal(detuning + dcShift);
            }

            /// 2. Populate the Off-Diagonals (Nearest-neighbor couplings)
            if constexpr (LevelCount > 1)
            {
                for (natural_t i = 0; i < LevelCount - 1; ++i)
                {
                    /// Coupling: Ω = -1/2 · μ_{i, i+1} · ε(t)
                    /// We extract μ from the full matrix row i and column i+1
                    const complex_t coupling = m_DipoleMatrix[i][i + 1] * (-0.5 * m_DriveAmplitude);

                    m_hamiltonianMatrix[SUPERDIAGONAL][i] = coupling;
                    m_hamiltonianMatrix[SUBDIAGONAL][i + 1] = coupling.conj();
                }
            }
        }
    };
}