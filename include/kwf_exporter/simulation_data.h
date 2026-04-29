#pragma once
#include <string>
#include "atomic_units.h"
#include "laser/laser.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    template <hilbert_space_t FullHilbertSpace>
    struct SimulationView
    {
        // Simulation m_time in SI
        real_t m_time;

        // Full 2D state vector
        StateVector<FullHilbertSpace> m_psi2D;

        // Bloch vectors
        complex_t m_alpha;
        complex_t m_beta;

        // Laser information in SI
        real_t m_laser1Wavelength;
        real_t m_laser1Intensity;
        real_t m_laser2Wavelength;
        real_t m_laser2Intensity;

        // Custom caption
        std::string m_title;
    };
}
