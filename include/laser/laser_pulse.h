#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Encapsulates the physical parameters of a driving laser field.
    struct LaserPulse
    {
        /// @brief Angular frequency of the light field (in a.u.).
        real_t m_omega = 0.0;
        /// @brief Electric field amplitude ε₀.
        real_t m_amplitude = 0.0; 
        /// @brief Temporal phase offset φ in radians.
        real_t m_phases = 0.0;     
    };

    /// @brief Encapsulates a two-photon driving configuration.
    /// @details Pairs the physical laser pulses (Pump and Stokes) with the specific 
    /// atomic energy level from which the multi-photon transition chain originates.
    struct TwoPhotonDrive
    {
        /// @brief The physical manifold index acting as the ground state for this transition chain.
        /// @note For Raman/Hadamard gates, this is typically Logical0Level. 
        /// For Rydberg transitions, it is shifted to Logical1Level.
        natural_t m_groundLevelOffset;

        /// @brief The pump laser driving the transition: |offset⟩ ↔ |offset + 1⟩.
        LaserPulse m_pump;

        /// @brief The Stokes laser driving the transition: |offset + 1⟩ ↔ |offset + 2⟩.
        LaserPulse m_stokes;
    };
}
