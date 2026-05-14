#pragma once
#include "core_types.h"
#include "constexprmath/constexpr_core_functions.h"
#include "constexprmath/constexpr_trigon.h"

/// @file Gate Library/Traits
/// @details The ideal unitary matrices are not included for direct calculations
/// but they will be used for gate fidelity and diagnostic purposes

namespace KetCat
{
	/// @brief Enumeration of supported quantum gate types.
    enum class GateType
    {
        I,

        X,
        Y,
        Z,
        H,

        RX,
        RY,
        RZ,

        CZ,
        CX,
        CCX,
        SWAP,
    };

    std::string gateNameToString(GateType gate)
    {
        switch (gate)
        {
            case GateType::I: return "Identity";
            case GateType::X: return "Pauli-X";
            case GateType::Y: return "Pauli-Y";
            case GateType::Z: return "Pauli-Z";
            case GateType::H: return "Hadamard";
            case GateType::RX: return "Rx";
            case GateType::RY: return "Ry";
            case GateType::RZ: return "Rz";
            case GateType::CZ: return "CPhase";
            case GateType::CX: return "CNOT";
            case GateType::CCX: return "Toffoli";
            case GateType::SWAP: return "SWAP";
            default: return "UnknownGate";
        }
	}


	/// @brief Compile-time traits and properties for each quantum gate type.
    template<GateType Gate>
    struct GateTraits;

    /// @brief Mathematical constants used for unitary gate definitions.
    constexpr real_t sqrt2 = 1.41421356237309505;
    constexpr real_t inv_sqrt2 = 1.0 / sqrt2;

    //========================================================
    // Identity
    //========================================================

    template<>
    struct GateTraits<GateType::I>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary() noexcept
        {
            return
            { {
                { complex_t(1.0, 0.0), complex_t(0.0, 0.0) },
                { complex_t(0.0, 0.0), complex_t(1.0, 0.0) }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "Identity";
        }
    };

    //========================================================
    // Pauli X
    //========================================================

    template<>
    struct GateTraits<GateType::X>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary() noexcept
        {
            return
            { {
                { complex_t(0.0, 0.0), complex_t(1.0, 0.0) },
                { complex_t(1.0, 0.0), complex_t(0.0, 0.0) }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "Pauli-X";
        }
    };

    //========================================================
    // Pauli Y
    //========================================================

    template<>
    struct GateTraits<GateType::Y>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary() noexcept
        {
            return
            { {
                { complex_t(0.0, 0.0),  complex_t(0.0, -1.0) },
                { complex_t(0.0, 1.0),  complex_t(0.0, 0.0) }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "Pauli-Y";
        }
    };

    //========================================================
    // Pauli Z
    //========================================================

    template<>
    struct GateTraits<GateType::Z>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary() noexcept
        {
            return
            { {
                { complex_t(1.0, 0.0),  complex_t(0.0, 0.0) },
                { complex_t(0.0, 0.0),  complex_t(-1.0, 0.0) }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "Pauli-Z";
        }
    };

    //========================================================
    // Hadamard
    //========================================================

    template<>
    struct GateTraits<GateType::H>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary() noexcept
        {
            constexpr real_t invsqrt2 =
                0.7071067811865475244;

            return
            { {
                {
                    complex_t(invsqrt2, 0.0),
                    complex_t(invsqrt2, 0.0)
                },
                {
                    complex_t(invsqrt2, 0.0),
                    complex_t(-invsqrt2, 0.0)
                }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "Hadamard";
        }
    };

    //========================================================
    // Rotation Y
    //========================================================

    template<>
    struct GateTraits<GateType::RY>
    {
        static constexpr natural_t QubitCount = 1;

        static constexpr matrix_t<2> unitary(real_t theta) noexcept
        {
            const real_t halfTheta =  theta / 2.0;

            return
            { {
                {
                    complex_t(
                        ConstexprMath::cos(halfTheta),
                        0.0
                    ),
                    complex_t(
                        -ConstexprMath::sin(halfTheta),
                        0.0
                    )
                },
                {
                    complex_t(
                        ConstexprMath::sin(halfTheta),
                        0.0
                    ),
                    complex_t(
                        ConstexprMath::cos(halfTheta),
                        0.0
                    )
                }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "RY";
        }
    };

    //========================================================
    // CNOT
    //========================================================

    template<>
    struct GateTraits<GateType::CX>
    {
        static constexpr natural_t QubitCount = 2;

        static constexpr matrix_t<4> unitary() noexcept
        {
            return
            { {
                {
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0)
                }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "CX";
        }
    };

    //========================================================
    // SWAP
    //========================================================

    template<>
    struct GateTraits<GateType::SWAP>
    {
        static constexpr natural_t QubitCount = 2;

        static constexpr matrix_t<4> unitary() noexcept
        {
            return
            { {
                {
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0)
                },
                {
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(0.0, 0.0),
                    complex_t(1.0, 0.0)
                }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "SWAP";
        }
    };

    //========================================================
    // Toffoli
    //========================================================

    template<>
    struct GateTraits<GateType::CCX>
    {
        static constexpr natural_t QubitCount = 3;

        static constexpr matrix_t<8>
            unitary() noexcept
        {
            return { {
                { complex_t(1,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(1,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(1,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(1,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(1,0), complex_t(0,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(1,0), complex_t(0,0), complex_t(0,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(1,0) },
                { complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(0,0), complex_t(1,0), complex_t(0,0) }
            } };
        }

        static constexpr const char* name() noexcept
        {
            return "CCX";
        }
    };
} 
