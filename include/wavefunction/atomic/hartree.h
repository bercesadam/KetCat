#pragma once
#include "constants.h"
#include "atomic_physics_core/elements.h"
#include "atomic_physics_core/quantum_number.h"
#include "atomic_physics_core/rydberg_quantum_defect.h"


namespace KetCat
{
    /// @brief Calculates the energy eigenvalue (Hartree energy) for a hydrogenic orbital.
    ///
    /// This computes the binding energy of the state, accounting for the 
    /// quantum defect effects in alkali-metal atoms. The energy is given 
    /// in atomic units (Hartree).
    ///
    /// The formula utilizes the Rydberg formula generalized with the effective 
    /// principal quantum number n* = n - δ_l:
    ///
    ///     E_n = - Z_eff² / (2 · (n - δ_l)²) = - Z_eff² / (2 · n*²)
    ///
    /// For Rydberg states of alkali metals, the effective nuclear charge Z_eff 
    /// is taken as 1, as the inner shell electrons screen the nucleus. For 
    /// pure hydrogenic systems, the quantum defect δ_l is zero, and the 
    /// formula reverts to the standard Bohr model energy levels.
    ///
    /// @param QNumbers 
    ///   QuantumNumber object containing the principal (n) and orbital (l) 
    ///   quantum numbers.
    /// @return 
    ///   The energy eigenvalue in Hartree units. Returns 0.0 if the effective 
    ///   principal quantum number is non-physical.
    template <quantum_number_t QuantumNumberType>
    constexpr real_t calculateHartreeEnergy(Element e, QuantumNumberType q) noexcept
    {
        // Use non-constexpr locals because these values depend on function parameters
        const real_t QuantumDefect = RydbergQuantumDefect::value(e, q);
        const real_t N_star = static_cast<real_t>(q.n()) - QuantumDefect;

        if (N_star <= 0.0)
        {
            return 0.0;
        }

		// For alkali metals, the effective nuclear charge Z_eff is approximately 1 due to screening by inner electrons.
        const real_t Z_eff = 1.0;
        return -(Z_eff * Z_eff) / (2.0 * N_star * N_star);
    }
}
