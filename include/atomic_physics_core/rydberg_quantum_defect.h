#pragma once
#include <utility>
#include "elements.h"
#include "quantum_number.h"


namespace KetCat
{
    /// @brief Quantum defects for alkali atoms (s, p, d, f states) based on experimental data.
    /// 
    /// @details The quantum defect δₗ modifies the effective principal quantum number n* = n - δₗ,
    /// Based on the NIST Atomic Spectra Database and other experimental sources, these quantum defects are essential for
    /// modeling the energy levels and wavefunctions of Rydberg states in alkali atoms, particularly for low-l states where the
    /// electron penetrates closer to the nucleus. This allows us to use the hydrogen wavefunction seeds with adjusted n*
    /// to approximate the Rydberg states of these atoms with reasonable accuracy,
    /// but it's more accurate for high-n states where the quantum defect becomes less significant.
    /// 
    /// Indexing: QuantumDefects[atom][l] gives the quantum defect for a given atom and orbital angular momentum l.
    class RydbergQuantumDefect
    {
        static constexpr std::array<std::array<real_t, 4>, 6> m_QuantumDefects =
        { {
                // s,     p,     d,     f
                {0.00,  0.00,  0.00,  0.00}, // H
                {0.40,  0.04,  0.00,  0.00}, // Li
                {1.35,  0.86,  0.01,  0.00}, // Na
                {2.18,  1.71,  0.28,  0.01}, // K
                {3.13,  2.65,  1.35,  0.02}, // Rb
                {4.05,  3.59,  2.47,  0.04}  // Cs
        }};

    public:
		/// @brief Get the quantum defect for a given atom, orbital angular momentum l, and principal quantum number n.
        template <quantum_number_t QuantumNumberType>
		static constexpr real_t value(Element element, QuantumNumberType q) noexcept
        {
            constexpr natural_t l = q.l();

            // For l >= 4, quantum defect is negligible
            if (l >= 4)
            {
				return 0.0; 
            }

            return m_QuantumDefects[std::to_underlying(element)][l];
		}
    };
}
