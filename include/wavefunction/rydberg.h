#pragma once
#include "quantum_number.h"

namespace KetCat
{
	///@brief Enumeration of alkali atoms
    enum class Atom : natural_t
    {
        H,
        Li,
        Na,
        K,
        Rb,
        Cs
    };

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
    class QuantumDefects
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
		constexpr real_t getEffectivePrincipalQuantumNumber(Atom atom, natural_t l, natural_t n) const noexcept
        {
            if (l >= 4)
            {
                // For l >= 4, quantum defect is negligible
				return n; 
            }

            real_t Delta = m_QuantumDefects[static_cast<natural_t>(atom)][l];
            return n - Delta;
		}
    };
}
