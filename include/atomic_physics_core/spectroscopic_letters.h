#pragma once
#include <concepts>
#include <type_traits>
#include "core_types.h"


namespace KetCat
{
	/// @brief Represents the azimuthal quantum number l, allowing us to use safe named constants for
    //  spectroscopic letters instead of raw integers.
    template<natural_t>
    struct AzimuthalQuantumNumber;

    template<> struct AzimuthalQuantumNumber<0> : std::integral_constant<natural_t, 0> {};
    template<> struct AzimuthalQuantumNumber<1> : std::integral_constant<natural_t, 1> {};
    template<> struct AzimuthalQuantumNumber<2> : std::integral_constant<natural_t, 2> {};
    template<> struct AzimuthalQuantumNumber<3> : std::integral_constant<natural_t, 3> {};
    template<> struct AzimuthalQuantumNumber<4> : std::integral_constant<natural_t, 4> {};
    template<> struct AzimuthalQuantumNumber<5> : std::integral_constant<natural_t, 5> {};
    template<> struct AzimuthalQuantumNumber<6> : std::integral_constant<natural_t, 6> {};
    template<> struct AzimuthalQuantumNumber<7> : std::integral_constant<natural_t, 7> {};
    template<> struct AzimuthalQuantumNumber<8> : std::integral_constant<natural_t, 8> {};
    template<> struct AzimuthalQuantumNumber<9> : std::integral_constant<natural_t, 9> {};

    /// @brief Type definitions for azimuthal quantum numbers corresponding to spectroscopic letters.
    /// Enclosing them in a namespace to not to pollute the main KetCat namespace with single-letter constants
    /// @group SpectroscopicLetters
    /// {
    namespace SpectroscopicLetters
    {
        using s = AzimuthalQuantumNumber<0>;
        using p = AzimuthalQuantumNumber<1>;
        using d = AzimuthalQuantumNumber<2>;
        using f = AzimuthalQuantumNumber<3>;
        using g = AzimuthalQuantumNumber<4>;
        using h = AzimuthalQuantumNumber<5>;
        using i = AzimuthalQuantumNumber<6>;
        using k = AzimuthalQuantumNumber<7>;
        using l = AzimuthalQuantumNumber<8>;
        using m = AzimuthalQuantumNumber<9>;
    }
    /// }

    /// @brief Type trait to check if a type is an AzimuthalQuantumNumber specialization
    template<typename T>
    struct is_azimuthal_quantum_number : std::false_type {};

    template<natural_t N>
    struct is_azimuthal_quantum_number<AzimuthalQuantumNumber<N>> : std::true_type {};
    
    /// @brief Concept to constrain template parameters to be valid azimuthal quantum number types
    template<typename T>
    concept azimuthal_quantum_number_t = is_azimuthal_quantum_number<T>::value;

    /// @brief Helper function to convert azimuthal quantum letter type
    /// to the corresponding spectroscopic letter for visualization reasons
    /// @param l AzimuthalQuantumNumber value
    /// @return The spectroscopic letter in char
    constexpr char getSpectroscopicLetter(natural_t l)
    {
        switch (l)
        {
            case 0: return 's';
            case 1: return 'p';
            case 2: return 'd';
            case 3: return 'f';
            case 4: return 'g';
            default: return '?';
        }
    }
}
