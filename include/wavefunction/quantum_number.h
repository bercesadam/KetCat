#pragma once
#include <concepts>
#include "core_types.h"

namespace KetCat
{
	/// @brief Represents the azimuthal quantum number l
    struct AzimuthalQuantumNumber
    {
		/// @brief The underlying value of the azimuthal quantum number l
        natural_t m_l;

		///@brief Constructor for AzimuthalQuantumNumber
        consteval AzimuthalQuantumNumber(natural_t l) : m_l(l){}
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

    /// @brief Representation of  quantum numbers (n, l, m)
    /// @details
    /// Encodes the triplet (n, l, m) with the usual constraints:
    ///  • n: principal quantum number, n ≥ 1 (controls energy and radial extent)  
    ///  • l: orbital angular momentum, 0 ≤ l ≤ n−1 (controls angular structure and parity)  
    ///  • m: magnetic quantum number, −l ≤ m ≤ l (controls angular orientation)  
    template<natural_t N, AzimuthalQuantumNumber L, int M = 0>
    class QuantumNumber
    {
        /// Compile-time checks for valid quantum numbers, replaces the old approach
        /// where I had a vast amount of predefined static constructors for each valid combination of (n, l, m)
        static_assert(N >= 1,
            "The principal quantum number (n) must be at least 1.");
        static_assert(L.m_l < N,
            "The azimuthal quantum number (l) must be less than the principal quantum number (n).");
        static_assert(M >= -static_cast<int>(L.m_l) && M <= static_cast<int>(L.m_l),
            "The magnetic quantum number (m) must satisfy -l ≤ m ≤ l.");

        /// @brief The underlying value of the principal quantum number n
        natural_t m_n;

        /// @brief The underlying value of the azimuthal quantum number l
        natural_t m_l;

        /// @brief The underlying value of the magnmetic quantum number m
        natural_t m_m;

    public:
		constexpr QuantumNumber() = default;

        ///@brief Getter for the principal quantum number n
        static constexpr natural_t n() noexcept { return N; }

        ///@brief Getter for the azimuthal quantum number l
        static constexpr natural_t l() noexcept { return L.m_l; }

        ///@brief Getter for the magetic quantum number m
        static constexpr int m() noexcept { return M; }

		/// @brief Calculate the energy eigenvalue in Hartree units for a hydrogenic orbital with quantum numbers (n, l, m).
        static constexpr real_t hartreeEnergy() noexcept {
            return -1.0 / (2.0 * N * N);
        }
    };

	/// @brief Concept to check if a type T has the required interface of quantum numbers (n, l, m) with appropriate return types.
    template<typename T>
    concept arithmetic_t = std::is_arithmetic_v<T>;

    template <typename T>
    concept quantum_number_t = requires
    {
        { T::n() } -> arithmetic_t;
        { T::l() } -> std::convertible_to<natural_t>;
        { T::m() } -> std::convertible_to<int>;
    };
}
