#pragma once
#include "core_types.h"

namespace KetCat
{
	/// @brief Slater's effective principal quantum numbers (n*)
	/// to be used for STO (Slater-type orbital) seed wavefunctions.
	class SlaterEffectiveQuantumNumber
    {
        ///@brief Slater's effective principal quantum numbers (n*) for alkali atoms, used to approximate the energy levels of Rydberg states.
        static constexpr std::array<real_t, 7> m_EffectivePrincipalQuantumNumbers =
        {
            0.0, // unused
            1.0, // n = 1
            2.0, // n = 2
            3.0, // n = 3
            3.7, // n = 4
            4.0, // n = 5
            4.2  // n = 6
        };

    public:
		/// @brief Get the effective principal quantum number n* for a given atom and principal quantum number n.
        constexpr real_t operator()(Atom atom, natural_t n) const noexcept
        {
            if (n >= 7)
            {
                // For n ≥ 7, the effective principal quantum number is approximately n - 2.5,
				// would correspond to Francium's (Fr) ground state, but on these higher energy levels we will
                // use Hidrogenic radials with Rydberg quantum defects.
				return n - 2.5;
            }

            return m_EffectivePrincipalQuantumNumbers[natural_t];
        }
    };
}
