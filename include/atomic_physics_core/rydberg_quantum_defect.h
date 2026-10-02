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
        // Quantum defect primary Ritz parameters (delta_0) for zero-field/low-l states.
        // Dimensions expanded to 11 to fit all 11 entries in the Element enum.
		// BIG WARNING: The quantum defects for alkali earth metals are for singly ionized states (X+), not neutral atoms,
        // added to support trapped ion simulations as well.
        static constexpr std::array<std::array<real_t, 4>, 11> m_QuantumDefects =
        { {
                // s,     p,     d,     f
                // --- Neutral Alkali Metals & Hydrogen ---
                // Ref: Drake, G. W. F. (2006). Springer Handbook of Atomic, Molecular, and Optical Physics. Springer.
                {0.00,  0.00,  0.00,  0.00}, // H  [Ref: Exact hydrogenic baseline, zero defect]
                {0.40,  0.04,  0.00,  0.00}, // Li [Ref: Lorenzen, C.-J., & Niemax, K. (1983). Phys. Scr., 27(4), 300]
                {1.35,  0.86,  0.01,  0.00}, // Na [Ref: Martin, W. C. (1980). J. Opt. Soc. Am., 70(1), 78-80]
                {2.18,  1.71,  0.28,  0.01}, // K  [Ref: Lorenzen, C.-J., & Niemax, K. (1983). Phys. Scr., 27(4), 300]
                {3.13,  2.65,  1.35,  0.02}, // Rb [Ref: Li, W., et al. (2003). Phys. Rev. A, 67(5), 052502]
                {4.05,  3.59,  2.47,  0.04}, // Cs [Ref: Weber, K.-H., & Sansonetti, C. J. (1987). Phys. Rev. A, 35(11), 4650]

                // --- Singly Ionized Alkaline-Earth Metals (X+) ---
                // Ref: Kramida, A., Ralchenko, Y., Reader, J., & NIST ASD Team (2024). NIST Atomic Spectra Database.
                {0.73,  0.36,  0.00,  0.00}, // Be+ [Ref: Johansson, L. (1962). Ark. Fys., 20, 489-498]
                {1.52,  1.10,  0.26,  0.01}, // Mg+ [Ref: Risberg, G. (1964). Ark. Fys., 28, 381-401]
                {2.34,  1.88,  0.89,  0.03}, // Ca+ [Ref: Edlén, B. (1956). Opt. Pura Apl., 10, 123; NIST ASD]
                {3.26,  2.81,  1.82,  0.07}, // Sr+ [Ref: Sansonetti, J. E., & Martin, W. C. (2005). J. Phys. Chem. Ref. Data, 34(4)]
                {4.15,  3.70,  2.71,  0.13}  // Ba+ [Ref: Curry, J. J. (2004). J. Phys. Chem. Ref. Data, 33(3), 725-745]
        } };

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
