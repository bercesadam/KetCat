#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "atomic_physics_core/atom.h"

namespace KetCat
{
    /// @brief Estimate the Slater orbital exponent (zeta) for a target valence shell
    ///        using Slater's shielding rules.
    ///
    /// This function approximates the effective nuclear charge felt by an electron
    /// in the selected shell/subshell, then converts it into the Slater-type orbital
    /// decay parameter:
    ///
    ///     zeta = Z_eff / n*
    ///
    /// where:
    ///   - Z_eff is the effective nuclear charge after electron shielding
    ///   - n*    is the effective principal quantum number (quantum defect corrected)
    ///
    /// This approximation is useful for:
    ///   - generating STO seed wavefunctions
    ///   - initializing eigenvalue solvers
    ///   - rough alkali atom orbital models
    ///
    /// @param elm               Chemical element (enum)
    /// @param targetShellIndex  Index of the shell/subshell in the Aufbau configuration
    ///                          for which zeta should be estimated
    /// @param nStar             Effective principal quantum number:
    ///                          n* = n - quantum_defect
    ///
    /// @return Estimated Slater exponent zeta
    constexpr real_t calculateZeta(Element elm, real_t nStar)
    {
		constexpr auto Config = Atom<elm>::getElectronConfiguration();

        // Atomic number = total number of electrons in a neutral atom.
        constexpr natural_t Z = std::to_underlying(elm);

        // Select the target subshell for which we estimate shielding.
        const ElectronShell Target = Config[targetShellIndex];

        // Total shielding contribution from all other electrons.
        real_t Shielding = 0.0;

        // Loop over every occupied shell/subshell in the atom.
        for (natural_t i = 0; i < MAX_SHELLS; ++i)
        {
            const ElectronShell Shell = Config[i];

            // Ignore unused entries in the fixed-size configuration array.
            if (Shell.m_numElectrons == 0)
            {
                continue;
            }

            // ------------------------------------------------------------
            // Slater shielding rules (simplified for ns / np valence shells)
            // ------------------------------------------------------------

            // Case 1:
            // Same subshell as the target electron.
            //
            // Other electrons in the same shell partially shield the nucleus.
            // According to Slater's rules:
            //   same-shell electrons contribute ~0.35 each
            //
            // Subtract 1 because the target electron does not shield itself.
            if (i == targetShellIndex)
            {
                Shielding += (Shell.m_numElectrons - 1) * 0.35;
            }

            // Case 2:
            // Electrons in the shell directly below the target shell (n - 1).
            //
            // These contribute stronger shielding:
            //   0.85 per electron
            else if (Shell.m_n == Target.m_n - 1)
            {
                Shielding += Shell.m_numElectrons * 0.85;
            }

            // Case 3:
            // Electrons in deeper inner shells (n - 2 and below).
            //
            // These almost completely screen the nucleus:
            //   1.00 per electron
            else if (Shell.m_n < Target.m_n - 1)
            {
                Shielding += Shell.m_numElectrons;
            }
        }

        // Effective nuclear charge experienced by the valence electron:
        //
        //     Z_eff = Z - shielding
        //
        const real_t zeff = static_cast<real_t>(Z) - shielding;

        // Slater orbital exponent:
        //
        //     zeta = Z_eff / n*
        //
        // Higher zeta:
        //   -> orbital is more compact / tightly bound
        //
        // Lower zeta:
        //   -> orbital is more diffuse
        return zeff / nStar;
    }

    template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
    struct SlaterOrbital
    {
        /// @param q      Kvantumszámok (itt n_star-t kellene használnod a sima n helyett)
        /// @param zeta   Az árnyékolt magtöltés paramétere (zeta = Z_eff / n_star)
        template <quantum_number_t QuantumNumberType>
        constexpr StateVector<HilbertSpace>
            operator()(QuantumNumberType q, Element elm) const noexcept
        {
            const double n_star = q.n_star();
            const natural_t l = q.l();

            StateVector<HilbertSpace> Psi{ complex_t::zero() };

            for (natural_t i = 1; i < HilbertSpace::Dim; ++i)
            {
                const double r = i * HilbertSpace::dx;

                // STO definíció: r^(n-1) * exp(-zeta * r)
                // De mivel te u(r) = r * R(r) formátumot használsz (reduced radial):
                // u(r) = r^(n_star) * exp(-zeta * r)

                const double rPow = std::pow(r, n_star);
                const double Exponential = ConstexprMath::exp<30>(-zeta * r);

                // Nincs Laguerre-polinom! Ez az STO lényege.
                const double Value = rPow * Exponential;
                Psi[i] = complex_t::fromReal(Value);
            }

            // A normalize() elintézi a bázis-független normalizációt.
            Psi.normalize();

            return Psi;
        }
    };
}