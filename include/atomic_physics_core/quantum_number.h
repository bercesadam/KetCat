#pragma once
#include <concepts>
#include <string>
#include "core_types.h"
#include "spectroscopic_letters.h"

namespace KetCat
{
    /// @brief Representation of  quantum numbers (n, l, m)
    /// @details
    /// Encodes the triplet (n, l, m) of quantum numbers as a compile-time type:
    ///  • n: principal quantum number, n ≥ 1 (controls energy and radial extent)  
    ///  • l: orbital angular momentum, 0 ≤ l ≤ n−1 (controls angular structure and parity)  
    ///  • m: magnetic quantum number, −l ≤ m ≤ l (controls angular orientation)  
    template<natural_t N, azimuthal_quantum_number_t L, int M = 0>
    class QuantumNumber
    {
        /// Compile-time checks for valid quantum numbers, replaces the old approach
        /// where I had a vast amount of predefined static constructors for each valid combination of (n, l, m)
        static_assert(N >= 1,
            "The principal quantum number (n) must be at least 1.");
        static_assert(L::value < N,
            "The azimuthal quantum number (l) must be less than the principal quantum number (n).");
        static_assert(M >= -static_cast<int>(L::value) && M <= static_cast<int>(L::value),
            "The magnetic quantum number (m) must satisfy -l ≤ m ≤ l.");

    public:
		constexpr QuantumNumber() = default;

        ///@brief Getter for the principal quantum number n
        static constexpr natural_t n() noexcept { return N; }

        ///@brief Getter for the azimuthal quantum number l
        static constexpr natural_t l() noexcept { return L::value; }

        ///@brief Getter for the magetic quantum number m
        static constexpr int m() noexcept { return M; }

		/// @brief Calculate the energy eigenvalue in Hartree units for a hydrogenic orbital with quantum numbers (n, l, m).
        static constexpr real_t hartreeEnergy() noexcept {
            return -1.0 / (2.0 * N * N);
        }

        static constexpr real_t rydbergEnergy() noexcept {
            return -1.0 / (2.0 * N * N);
		}
    };

	/// @brief Concept to check if a type T has the required interface of quantum numbers (n, l, m) with appropriate return types.
    template<typename T>
    concept arithmetic_t = std::is_arithmetic_v<T>;

    template <typename T>
    concept quantum_number_t = requires
    {
        { T::n() } -> std::convertible_to<natural_t>;
        { T::l() } -> std::convertible_to<natural_t>;
        { T::m() } -> std::convertible_to<int>;
    };

    template <quantum_number_t QNumber>
    std::string quantumNumberToString()
    {
        return std::to_string(QNumber::n()) + getSpectroscopicLetter(QNumber::l());
    }
}
