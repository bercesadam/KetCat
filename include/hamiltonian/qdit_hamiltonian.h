#pragma once
#include "core_types.h"


namespace KetCat
{
    enum class RotationAxis
    {
        X,
        Y,
        Z
    };

    /// @brief Driven qudit Hamiltonian on a reduced energy space.
    /// 
    /// - Free H: H_ii = E_i  (actual eigenenergies in Hartree, e.g. from seeds)
    /// - X/Y drive: single-tone at ω_d = E_b - E_a
    /// - Leakage ladder: weak global nearest-neighbor coupling, driven at ω_d
    ///
    /// NOTE: tridiagonal representation only supports nearest-neighbor coupling.
    /// Therefore, (a,b) must be adjacent: b = a+1.
    template<natural_t Dim>
    class SingleQditGateHamiltonian
    {
        // Energies of the reduced basis (Hartree)
        std::array<real_t, Dim> m_E;

        // Drive envelope parameters
        real_t m_Omega;
        real_t m_lambda;
        real_t m_leakageFactor;

        // Logical pair (must be adjacent for tridiagonal coupling)
        natural_t m_a;
        natural_t m_b;

        RotationAxis m_axis;

    public:
        /// @param axis           X, Y, or Z rotation frame (Z = free Hamiltonian only)
        /// @param energies       Physical energies of the reduced basis (Hartree)
        /// @param Omega          Rabi-rate scale for the logical drive (a.u.)
        /// @param lambda         Characteristic leakage-decay scale across ladder
        /// @param leakageFactor  Global scaling of leakage vs logical drive
        /// @param logicalA       Lower logical level index (0-based)
        /// @param logicalB       Upper logical level index (0-based), must be A+1
        constexpr SingleQditGateHamiltonian(
            RotationAxis axis,
            const std::array<real_t, Dim>& energies,
            real_t Omega,
            real_t lambda,
            real_t leakageFactor,
            natural_t logicalA,
            natural_t logicalB
        ) noexcept
            : m_E(energies)
            , m_Omega(Omega)
            , m_lambda(lambda)
            , m_leakageFactor(leakageFactor)
            , m_a(logicalA < logicalB ? logicalA : logicalB)
            , m_b(logicalA < logicalB ? logicalB : logicalA)
            , m_axis(axis)
        {}

        constexpr tridiagonal_matrix_t<Dim>
        operator()(real_t t) const noexcept
        {
            tridiagonal_matrix_t<Dim> H{};

            /// ====================================================
            /// Free Hamiltonian: use actual energies
            /// ====================================================
            for (natural_t n = 0; n < Dim; ++n)
            {
                H[MAINDIAGONAL][n] = complex_t::fromReal(m_E[n]);
            }

            if (m_axis == RotationAxis::Z)
            {
                // Z "rotation": just evolve under free Hamiltonian
                return H;
            }

            /// ====================================================
            /// Resonant logical drive (single tone)
            /// ====================================================
            // Require nearest-neighbor logical pair for tridiagonal coupling
            if (!(m_b == m_a + 1))
            {
                // In production, consider throwing or guarding at construction.
                // Here, we just early-return the free H if invalid pair.
                // assert(false && "XYQditGateHamiltonian: logical pair must be adjacent (b == a+1)");
                return H;
            }

            const real_t omega_d = m_E[m_b] - m_E[m_a];      // Hartree (ħ = 1)
            const real_t drive   = m_Omega * ConstexprMath::cos(omega_d * t);

            if (m_axis == RotationAxis::X)
            {
                // <a|H|b> = <b|H|a> = Ω cos(ω_d t)
                H[SUPERDIAGONAL][m_a] += complex_t::fromReal(drive);     // (a, a+1)
                H[SUBDIAGONAL][m_b]   += complex_t::fromReal(drive);     // (b, b-1)
            }
            else if (m_axis == RotationAxis::Y)
            {
                // <a|H|b> = -i Ω cos(ω_d t),  <b|H|a> = +i Ω cos(ω_d t)
                H[SUPERDIAGONAL][m_a] += complex_t(0.0, -drive);
                H[SUBDIAGONAL][m_b]   += complex_t(0.0,  drive);
            }

            /// ====================================================
            /// Weak global leakage ladder drive (nearest-neighbor)
            /// ====================================================
            for (natural_t n = 0; n + 1 < Dim; ++n)
            {
                if ((n == m_a) && (n + 1 == m_b))
                    continue;

                // Envelope: scaled HO-like √(n+1) times exponential falloff
                const real_t coupling =
                    m_leakageFactor *
                    m_Omega *
                    ConstexprMath::sqrt(real_t(n + 1)) *
                    ConstexprMath::exp<10>(-real_t(n) / m_lambda);

                // Single-tone global drive at ω_d:
                // Off-resonant couplings are naturally suppressed by detuning
                // via the free Hamiltonian energies on the diagonal.
                const real_t leakageDrive = coupling * ConstexprMath::cos(omega_d * t);

                if (m_axis == RotationAxis::X)
                {
                    H[SUPERDIAGONAL][n]     += complex_t::fromReal(leakageDrive);
                    H[SUBDIAGONAL][n + 1]   += complex_t::fromReal(leakageDrive);
                }
                else if (m_axis == RotationAxis::Y)
                {
                    H[SUPERDIAGONAL][n]     += complex_t(0.0, -leakageDrive);
                    H[SUBDIAGONAL][n + 1]   += complex_t(0.0,  leakageDrive);
                }
            }

            return H;
        }
    };
}
