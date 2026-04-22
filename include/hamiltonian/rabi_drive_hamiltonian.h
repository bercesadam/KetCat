#pragma once 
#include "core_types.h"

namespace KetCat
{

    
    real_t safeDetuning(real_t delta) noexcept
    {
        constexpr real_t eps = 1e-10;
        return (ConstexprMath::abs(delta) > eps) ? delta : (delta >= 0 ? eps : -eps);
    }


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

            const real_t E2 = m_DriveAmplitude * m_DriveAmplitude;

            // ------------------------------------------------------------
            // 1. Main diagonal: rotating-frame detuning + AC Stark shift
            // ------------------------------------------------------------
            for (natural_t i = 0; i < LevelCount; ++i)
            {
                // Bare rotating-frame detuning
                real_t detuning = m_Energies[i] - static_cast<real_t>(i) * m_DriveOmega;

                // AC Stark shift (second order, nearest neighbours only)
                real_t acStark = 0.0;

                // Coupling to upper neighbour (i -> i+1)
                if (i + 1 < LevelCount)
                {
                    const real_t delta =
                        (m_Energies[i + 1] - m_Energies[i]) - m_DriveOmega;

                    const real_t OmegaSq =
                        m_DipoleMatrix[i][i + 1].normSquared() * E2;

                    acStark -= OmegaSq / (4.0 * safeDetuning(delta));
                }

                // Coupling to lower neighbour (i -> i-1)
                if (i > 0)
                {
                    const real_t delta =
                        (m_Energies[i - 1] - m_Energies[i]) + m_DriveOmega;

                    const real_t OmegaSq =
                        m_DipoleMatrix[i][i - 1].normSquared() * E2;

                    acStark -= OmegaSq / (4.0 * safeDetuning(delta));
                }

                m_hamiltonianMatrix[MAINDIAGONAL][i] =
                    complex_t::fromReal(detuning);// + acStark);
            }

            // ------------------------------------------------------------
            // 2. Off-diagonals: Rabi couplings (unchanged)
            // ------------------------------------------------------------
            if constexpr (LevelCount > 1)
            {
                for (natural_t i = 0; i < LevelCount - 1; ++i)
                {
                    const complex_t coupling =
                        m_DipoleMatrix[i][i + 1]
                        * (-0.5 * m_DriveAmplitude);

                    m_hamiltonianMatrix[SUPERDIAGONAL][i] = coupling;
                    m_hamiltonianMatrix[SUBDIAGONAL][i + 1] = coupling.conj();
                }
            }
        }
    };
}