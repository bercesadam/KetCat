#pragma once
#include "core_types.h"

namespace KetCat
{
    /// @brief Compact representation of hydrogenic quantum numbers (n, l, m) in spectroscopic notation.
    /// @details
    /// Encodes the triplet (n, l, m) with the usual constraints:
    ///  • n: principal quantum number, n ≥ 1 (controls energy and radial extent)  
    ///  • l: orbital angular momentum, 0 ≤ l ≤ n−1 (controls angular structure and parity)  
    ///  • m: magnetic quantum number, −l ≤ m ≤ l (controls angular orientation)  
    ///
    /// Spectroscopic letters: s≡l=0, p≡1, d≡2, f≡3, g≡4, …  
    /// This type offers a set of `constexpr` named constructors for valid (n,l,m) choices to
    /// avoid selecting invalid configurations.
    class QuantumNumber
    {
        unsigned int m_n; // principal quantum number
        unsigned int m_l; // orbital angular momentum
        int m_m;          // magnetic quantum number

        /// @brief Construct a quantum number triple (n,l,m).
        /// @param n Principal quantum number (n ≥ 1).
        /// @param l Orbital angular momentum (0 ≤ l ≤ n−1).
        /// @param m Magnetic quantum number (−l ≤ m ≤ l).
        constexpr QuantumNumber(unsigned int n, unsigned int l, int m)
            : m_n(n), m_l(l), m_m(m)
        {
        }

    public:
        /// @brief Principal quantum number n.
        constexpr unsigned int n() const { return m_n; }

        /// @brief Orbital angular momentum l.
        constexpr unsigned int l() const { return m_l; }

        /// @brief Magnetic quantum number m.
        constexpr unsigned int m() const { return m_m; }

        /// @brief Hydrogenic energy in Hartree units for the given n.
        /// @details Eₙ = −1 / (2 n²). Assumes non‑relativistic Coulomb problem (Z=1).
        constexpr real_t hartreeEnergy() const noexcept
        {
            return -1.0 / (2.0 * m_n * m_n);
        }

        // --------------------------
        // Static constexpr constructors for allowed (n,l,m)
        // --------------------------

        /// @name Constructors for (n,l) only with m = 0 (eg. to be used in 1D examples)
        /// @{
        static constexpr QuantumNumber _1s() { return QuantumNumber{ 1, 0, 0 }; }
        static constexpr QuantumNumber _2s() { return QuantumNumber{ 2, 0, 0 }; }
        static constexpr QuantumNumber _2p() { return QuantumNumber{ 2, 1, 0 }; }
        static constexpr QuantumNumber _3s() { return QuantumNumber{ 3, 0, 0 }; }
        static constexpr QuantumNumber _3p() { return QuantumNumber{ 3, 1, 0 }; }
        static constexpr QuantumNumber _3d() { return QuantumNumber{ 3, 2, 0 }; }
        static constexpr QuantumNumber _4s() { return QuantumNumber{ 4, 0, 0 }; }
        static constexpr QuantumNumber _4p() { return QuantumNumber{ 4, 1, 0 }; }
        static constexpr QuantumNumber _4d() { return QuantumNumber{ 4, 2, 0 }; }
        static constexpr QuantumNumber _4f() { return QuantumNumber{ 4, 3, 0 }; }
        static constexpr QuantumNumber _5g() { return QuantumNumber{ 5, 4, 0 }; }
        static constexpr QuantumNumber _6h() { return QuantumNumber{ 6, 5, 0 }; }
        static constexpr QuantumNumber _7i() { return QuantumNumber{ 7, 6, 0 }; }
        static constexpr QuantumNumber _8k() { return QuantumNumber{ 8, 7, 0 }; }
        /// @}

        /// @name n = 1
        /// @{
        static constexpr QuantumNumber _1s_m0() { return QuantumNumber{ 1, 0, 0 }; }
        /// @}

        /// @name n = 2 (s, p with m ∈ {−l,…,l})
        /// @{
        static constexpr QuantumNumber _2s_m0()      { return QuantumNumber{ 2, 0, 0 }; }
        static constexpr QuantumNumber _2p_m_minus1(){ return QuantumNumber{ 2, 1, -1 }; }
        static constexpr QuantumNumber _2p_m0()      { return QuantumNumber{ 2, 1,  0 }; }
        static constexpr QuantumNumber _2p_m1()      { return QuantumNumber{ 2, 1,  1 }; }
        /// @}

        /// @name n = 3 (s, p, d with m ∈ {−l,…,l})
        /// @{
        static constexpr QuantumNumber _3s_m0()       { return QuantumNumber{ 3, 0, 0 }; }
        static constexpr QuantumNumber _3p_m_minus1() { return QuantumNumber{ 3, 1, -1 }; }
        static constexpr QuantumNumber _3p_m0()       { return QuantumNumber{ 3, 1,  0 }; }
        static constexpr QuantumNumber _3p_m1()       { return QuantumNumber{ 3, 1,  1 }; }
        static constexpr QuantumNumber _3d_m_minus2() { return QuantumNumber{ 3, 2, -2 }; }
        static constexpr QuantumNumber _3d_m_minus1() { return QuantumNumber{ 3, 2, -1 }; }
        static constexpr QuantumNumber _3d_m0()       { return QuantumNumber{ 3, 2,  0 }; }
        static constexpr QuantumNumber _3d_m1()       { return QuantumNumber{ 3, 2,  1 }; }
        static constexpr QuantumNumber _3d_m2()       { return QuantumNumber{ 3, 2,  2 }; }
        /// @}

        /// @name n = 4 (s, p, d, f with m ∈ {−l,…,l})
        /// @{
        static constexpr QuantumNumber _4s_m0()       { return QuantumNumber{ 4, 0, 0 }; }
        static constexpr QuantumNumber _4p_m_minus1() { return QuantumNumber{ 4, 1, -1 }; }
        static constexpr QuantumNumber _4p_m0()       { return QuantumNumber{ 4, 1,  0 }; }
        static constexpr QuantumNumber _4p_m1()       { return QuantumNumber{ 4, 1,  1 }; }
        static constexpr QuantumNumber _4d_m_minus2() { return QuantumNumber{ 4, 2, -2 }; }
        static constexpr QuantumNumber _4d_m_minus1() { return QuantumNumber{ 4, 2, -1 }; }
        static constexpr QuantumNumber _4d_m0()       { return QuantumNumber{ 4, 2,  0 }; }
        static constexpr QuantumNumber _4d_m1()       { return QuantumNumber{ 4, 2,  1 }; }
        static constexpr QuantumNumber _4d_m2()       { return QuantumNumber{ 4, 2,  2 }; }
        static constexpr QuantumNumber _4f_m_minus3() { return QuantumNumber{ 4, 3, -3 }; }
        static constexpr QuantumNumber _4f_m_minus2() { return QuantumNumber{ 4, 3, -2 }; }
        static constexpr QuantumNumber _4f_m_minus1() { return QuantumNumber{ 4, 3, -1 }; }
        static constexpr QuantumNumber _4f_m0()       { return QuantumNumber{ 4, 3,  0 }; }
        static constexpr QuantumNumber _4f_m1()       { return QuantumNumber{ 4, 3,  1 }; }
        static constexpr QuantumNumber _4f_m2()       { return QuantumNumber{ 4, 3,  2 }; }
        static constexpr QuantumNumber _4f_m3()       { return QuantumNumber{ 4, 3,  3 }; }
        /// @}
    };
}
