#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "atomic_physics_core/atom.h"
#include "atomic_physics_core/slater_effective_n.h"


namespace KetCat
{
    /// @brief Helper function to estimate the Slater orbital exponent (zeta) for a
    ///        target valence shell using Slater's shielding rules.
    ///
    /// This function approximates the effective nuclear charge felt by an electron
    /// in the selected shell/subshell, then converts it into the Slater-type orbital
    /// decay parameter:
    ///
    ///     ζ = Z_eff / n*
    ///
    /// where:
    ///   - Z_eff is the effective nuclear charge after electron shielding
    ///   - n*    is the effective principal quantum number (according to Slater's rules)
    ///
    /// This approximation is useful for:
    ///   - generating STO seed wavefunctions
    ///   - initializing eigenvalue solvers
    ///   - rough alkali atom orbital models
    ///
    /// @param element               Chemical element (enum)
    /// @param targetShellIndex  Index of the shell/subshell in the Aufbau configuration
    ///                          for which zeta should be estimated
    /// @param nStar             Effective principal quantum number:
    ///                          n* = n - quantum_defect
    ///
    /// @return Estimated Slater exponent zeta
    template <Element element>
    constexpr real_t calculateZeta(real_t nStar)
    {
        constexpr natural_t TargetShellIndex = Atom<element>::getOuterShellIndex();
		constexpr auto Config = Atom<element>::getElectronConfiguration();

        // Atomic number = total number of electrons in a neutral atom.
		constexpr natural_t Z = AtomicNumber<element>::value;

        // Select the target subshell for which we estimate shielding.
        const ElectronShell Target = Config[TargetShellIndex];

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
            if (i == TargetShellIndex)
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
        const real_t Z_eff = static_cast<real_t>(Z) - Shielding;

        // Slater orbital exponent:
        //
        //     ζ = Z_eff / n*
        //
        // Higher zeta means that the orbital is more compact / tightly bound
        // Lower zeta means that the orbital is more diffuse
        return Z_eff / nStar;
    }

    
    /// @brief Construct a Slater-type reduced radial wavefunction seed u(r) flattened to 1D.
    ///
    /// @details
    ///   This implements a Slater-Type Orbital (STO) inspired radial seed in reduced-radial form.
    ///   The returned function is u(r) = r · R(r), discretized on a 1D radial grid and
    ///   normalized on that grid.
    ///
    ///   Unlike hydrogenic orbitals, this form avoids Laguerre polynomials and instead uses
    ///   a simple exponential decay with an effective screening parameter ζ. The construction
    ///   is designed to preserve the correct near-origin behavior required by the radial
    ///   Schrödinger equation while remaining numerically simple and stable.
    ///
    ///   The reduced radial function is defined as:
    ///
    ///     u(r) = r^{ℓ+1} · exp(−(ζ / n*) · r)
    ///
    ///   where:
    ///     - ℓ    is the orbital angular momentum quantum number
    ///     - ζ    is the screened nuclear charge parameter (from Slater rules)
    ///     - n*   is the effective principal quantum number (n corrected by quantum defect)
    ///
    ///   Key properties:
    ///
    ///     • Near-origin behavior:
    ///         u(r) ~ r^{ℓ+1}
    ///       This satisfies the boundary condition imposed by the centrifugal term
    ///       ℓ(ℓ+1)/r² in the radial Schrödinger equation.
    ///
    ///     • Radial extent scaling:
    ///         exp(−(ζ / n*) r)
    ///       The inclusion of n* in the decay term increases the spatial extent of
    ///       higher-n states, mimicking the physical scaling of orbital size.
    ///
    ///   The full spatial wavefunction is:
    ///
    ///     ψ_{nℓm}(r,θ,φ) = (u_{nℓ}(r) / r) · Y_{ℓm}(θ,φ)
    ///
    ///   where Y_{ℓm} are spherical harmonics. The Hamiltonian is m-independent;
    ///   the magnetic quantum number m enters only through the angular part.
    ///
    ///   This construction is particularly useful for:
    ///     - initializing radial eigenvalue solvers
    ///     - generating physically consistent seed states
    ///     - approximate modeling of alkali valence orbitals
    ///
    /// @tparam HilbertSpace
    ///   Discrete 1D spatial Hilbert space defining the radial grid size and spacing.
    template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace, Element element>
    struct SlaterOrbitalRadial
    {
        /// @brief Generate a Slater-type reduced radial orbital.
        ///
        /// @param element   Chemical element descriptor (e.g. Li, Na, K)
		/// @param q         Quantum number pair (n, l)
        ///
        /// @return
        ///   Normalized state vector representing u(r)
        template <quantum_number_t QuantumNumberType>
        constexpr StateVector<HilbertSpace>
            operator()(QuantumNumberType q) const noexcept
        {
            // Effective (possibly non-integer) principal quantum number n*
            const real_t l = static_cast<real_t>(q.l());
            const real_t N_star = SlaterEffectiveQuantumNumber::value(q.n());

            // Estimate the Slater orbital exponent ζ using effective nuclear charge and n*.
            const real_t Zeta = calculateZeta<element>(N_star);

            StateVector<HilbertSpace> Psi{ complex_t::zero() };

            for (natural_t i = 1; i < HilbertSpace::Dim; ++i)
            {
                const real_t r = HilbertSpace::gridToR(i);

                const real_t nFactor = N_star - (l + 1.0);
                const real_t logValue =
                    (l + 1.0 + 0.3 * nFactor) * log(r) - Zeta * r;

               // const real_t logValue = (l + 1.0) * ConstexprMath::log(r) - (Zeta / N_star) * r;
                Psi[i] = complex_t::fromReal(ConstexprMath::exp(logValue));
            }

            // Enforce discrete radial normalization: Σ |u|² · Δr = 1
            Psi.normalize();

            return Psi;
        }
    };
}
