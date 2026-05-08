#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Encapsulates the physical parameters of a driving laser field.
    struct LaserPulse
    {
        real_t m_omega = 0.0;     ///< Angular frequency of the light field (in a.u.).
        real_t m_amplitude = 0.0; ///< Electric field amplitude ε₀.
        real_t m_phase = 0.0;     ///< Temporal phase offset φ in radians.
    };
}
