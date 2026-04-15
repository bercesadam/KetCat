#pragma once
#include <type_traits>
#include "core_types.h"


namespace KetCat
{
	///@brief Enumeration to represent and identify alkali atoms
    enum class Element : natural_t
    {
        H,
        Li,
        Na,
        K,
        Rb,
        Cs,
    };

    /// @brief Tag struct to associate each Element with its corresponding atomic number Z.
    template<Element>
    struct AtomicNumber; // intentionally undefined

    /// @brief Specializations of the AtomicNumber struct for each Element, providing the atomic number as a compile-time constant.
     /// @group AtomicNumbers
    /// {
    template<> struct AtomicNumber<Element::H>  : std::integral_constant<natural_t, 1>  {};
    template<> struct AtomicNumber<Element::Li> : std::integral_constant<natural_t, 3>  {};
    template<> struct AtomicNumber<Element::Na> : std::integral_constant<natural_t, 11> {};
    template<> struct AtomicNumber<Element::K>  : std::integral_constant<natural_t, 19> {};
    template<> struct AtomicNumber<Element::Rb> : std::integral_constant<natural_t, 37> {};
    template<> struct AtomicNumber<Element::Cs> : std::integral_constant<natural_t, 55> {};
    /// }

    
    /// @brief Type trait to check if a type is an AtomicNumber specialization
    template<typename T>
    struct is_atomic_number : std::false_type {};

    template<Element E>
    struct is_atomic_number<AtomicNumber<E>> : std::true_type {};

    /// @brief Concept to constrain template parameters to be valid atomic number types (i.e., specializations of AtomicNumber)
    template<typename T>
    concept atomic_number_t = is_atomic_number<T>::value;
}
