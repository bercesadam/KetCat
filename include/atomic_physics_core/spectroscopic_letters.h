#pragma once
#include <concepts>
#include "core_types.h"

namespace KetCat
{
	/// @brief Represents the azimuthal quantum number l, allowing us to use safe named constants for
    //  spectroscopic letters instead of raw integers.
    struct AzimuthalQuantumNumber
    {
        /// @brief The underlying value of the azimuthal quantum number l
        natural_t m_l;

        ///@brief Constructor for AzimuthalQuantumNumber
        consteval AzimuthalQuantumNumber(natural_t l) : m_l(l) {}
    };

    /// @brief Named constants for azimuthal quantum numbers corresponding to spectroscopic letters.
    /// Enclosing them in a namespace to not to pollute the main KetCat namespace with single-letter constants
    /// @group SpectroscopicLetters
    /// {
    namespace SpectroscopicLetters
    {
        constexpr AzimuthalQuantumNumber s{ 0 };
        constexpr AzimuthalQuantumNumber p{ 1 };
        constexpr AzimuthalQuantumNumber d{ 2 };
        constexpr AzimuthalQuantumNumber f{ 3 };
        constexpr AzimuthalQuantumNumber g{ 4 };
        constexpr AzimuthalQuantumNumber h{ 5 };
        constexpr AzimuthalQuantumNumber i{ 5 };
        constexpr AzimuthalQuantumNumber k{ 6 };
        constexpr AzimuthalQuantumNumber l{ 7 };
        constexpr AzimuthalQuantumNumber m{ 8 };
        constexpr AzimuthalQuantumNumber n{ 9 };
    }
    /// }
}
