#pragma once
#include "core_types.h"

namespace KetCat
{
    ;
    ///@brief Represents a triplet of quantum numbers (n, l, m) using the spectroscopic notation:
    ///- n : Principal quantum number (n >= 1), controls energy and radial extent
    ///- l : Orbital angular momentum quantum number (0 <= l <= n-1), controls angular structure and parity
    ///- m : Magnetic quantum number (-l <= m <= l), controls angular orientation
    class QuantumNumber
    {
        unsigned int m_n; // principal quantum number
        unsigned int m_l; // orbital angular momentum
        int m_m; // magnetic quantum number

        constexpr QuantumNumber(unsigned int n, unsigned int l, int m)
            : m_n(n), m_l(l), m_m(m)
        {
        }

    public:
        constexpr unsigned int n() const { return m_n; }
        constexpr unsigned int l() const { return m_l; }
        constexpr unsigned int m() const { return m_m; }

        constexpr real_t hartreeEnergy() noexcept
        {
            return -1.0 / (2.0 * m_n * m_n);
        }

        // --------------------------
        // Static constexpr constructors for all allowed n,l,m
        // --------------------------

        // Legacy constructors for n,l only (m=0)
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

        // n=1
        static constexpr QuantumNumber _1s_m0() { return QuantumNumber{ 1, 0, 0 }; }

        // n=2
        static constexpr QuantumNumber _2s_m0() { return QuantumNumber{ 2, 0, 0 }; }
        static constexpr QuantumNumber _2p_m_minus1() { return QuantumNumber{ 2, 1,-1 }; }
        static constexpr QuantumNumber _2p_m0() { return QuantumNumber{ 2, 1, 0 }; }
        static constexpr QuantumNumber _2p_m1() { return QuantumNumber{ 2, 1, 1 }; }

        // n=3
        static constexpr QuantumNumber _3s_m0() { return QuantumNumber{ 3, 0, 0 }; }
        static constexpr QuantumNumber _3p_m_minus1() { return QuantumNumber{ 3, 1,-1 }; }
        static constexpr QuantumNumber _3p_m0() { return QuantumNumber{ 3, 1, 0 }; }
        static constexpr QuantumNumber _3p_m1() { return QuantumNumber{ 3, 1, 1 }; }
        static constexpr QuantumNumber _3d_m_minus2() { return QuantumNumber{ 3, 2,-2 }; }
        static constexpr QuantumNumber _3d_m_minus1() { return QuantumNumber{ 3, 2,-1 }; }
        static constexpr QuantumNumber _3d_m0() { return QuantumNumber{ 3, 2, 0 }; }
        static constexpr QuantumNumber _3d_m1() { return QuantumNumber{ 3, 2, 1 }; }
        static constexpr QuantumNumber _3d_m2() { return QuantumNumber{ 3, 2, 2 }; }

        // n=4
        static constexpr QuantumNumber _4s_m0() { return QuantumNumber{ 4, 0, 0 }; }
        static constexpr QuantumNumber _4p_m_minus1() { return QuantumNumber{ 4, 1,-1 }; }
        static constexpr QuantumNumber _4p_m0() { return QuantumNumber{ 4, 1, 0 }; }
        static constexpr QuantumNumber _4p_m1() { return QuantumNumber{ 4, 1, 1 }; }
        static constexpr QuantumNumber _4d_m_minus2() { return QuantumNumber{ 4, 2,-2 }; }
        static constexpr QuantumNumber _4d_m_minus1() { return QuantumNumber{ 4, 2,-1 }; }
        static constexpr QuantumNumber _4d_m0() { return QuantumNumber{ 4, 2, 0 }; }
        static constexpr QuantumNumber _4d_m1() { return QuantumNumber{ 4, 2, 1 }; }
        static constexpr QuantumNumber _4d_m2() { return QuantumNumber{ 4, 2, 2 }; }
        static constexpr QuantumNumber _4f_m_minus3() { return QuantumNumber{ 4, 3,-3 }; }
        static constexpr QuantumNumber _4f_m_minus2() { return QuantumNumber{ 4, 3,-2 }; }
        static constexpr QuantumNumber _4f_m_minus1() { return QuantumNumber{ 4, 3,-1 }; }
        static constexpr QuantumNumber _4f_m0() { return QuantumNumber{ 4, 3, 0 }; }
        static constexpr QuantumNumber _4f_m1() { return QuantumNumber{ 4, 3, 1 }; }
        static constexpr QuantumNumber _4f_m2() { return QuantumNumber{ 4, 3, 2 }; }
        static constexpr QuantumNumber _4f_m3() { return QuantumNumber{ 4, 3, 3 }; }
    };
}