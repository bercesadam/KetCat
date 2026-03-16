#pragma once 
#include "core_types.h"

namespace KetCat
{

    /// @brief Rotation axes for single-qudit gates.
    ///
    /// Defines along which axis the qudit rotation (drive) acts.
    enum class RotationAxis
    {
        X, ///< Rotation around X axis (real dipole coupling)
        Y, ///< Rotation around Y axis (imaginary dipole coupling)
        Z  ///< Rotation around Z axis (phase shift only)
    };

    /// @brief Single-qudit gate Hamiltonian for a drive between two levels.
    ///
    /// @tparam Dim Dimension of the qudit (number of energy levels)
    ///
    /// ------------------------------------------------------------------------
    /// Mathematical construction
    /// ------------------------------------------------------------------------
    ///
    /// The Hamiltonian in the energy eigenbasis is tridiagonal with elements:
    ///
    ///      H_{ii} = E_i  (diagonal energies)
    ///      H_{i,i+1} = Ω * cos(ω_d t)  (off-diagonal dipole coupling)
    ///
    /// for levels |a⟩ ↔ |b⟩ and drive frequency:
    ///
    ///      ω_d = (E_b - E_a) + Δ
    ///
    /// Depending on rotation axis:
    ///
    ///  - X: real coupling
    ///  - Y: imaginary coupling
    ///  - Z: phase shift on diagonal elements only
    ///
    /// ------------------------------------------------------------------------
    /// Physical meaning
    /// ------------------------------------------------------------------------
    ///
    /// This object models a coherent drive between two qudit levels,
    /// generating the time-dependent Hamiltonian for simulation.
    template<natural_t Dim>
    class SingleQuditGateHamiltonian
    {
    private:

        /// @brief Pointer to array of energy eigenvalues {E_0, ..., E_{Dim-1}}
        const std::array<real_t, Dim> m_Energies;

        /// @brief Lower level of the driven transition |a⟩
        natural_t m_a;

        /// @brief Upper level of the driven transition |b⟩
        natural_t m_b;

        /// @brief Drive amplitude Ω
        real_t m_Omega;

        /// @brief Detuning Δ from exact resonance
        real_t m_Detuning;

        /// @brief Rotation axis of the drive (X, Y, Z)
        RotationAxis m_axis;

    public:

        /// @brief Construct a single-qudit gate Hamiltonian.
        ///
        /// @param energies Energy eigenvalues of the qudit levels
        /// @param a        Lower level of the driven transition
        /// @param b        Upper level of the driven transition
        /// @param Omega    Drive amplitude
        /// @param detuning Detuning from resonance
        /// @param axis     Rotation axis (X, Y, Z)
        constexpr SingleQuditGateHamiltonian(
            const std::array<real_t, Dim> energies,
            natural_t a,
            natural_t b,
            real_t Omega,
            real_t detuning,
            RotationAxis axis) noexcept
            :
            m_Energies(energies),
            m_a(a),
            m_b(b),
            m_Omega(Omega),
            m_Detuning(detuning),
            m_axis(axis)
        {}

        /// @brief Generate the time-dependent Hamiltonian matrix at time t.
        ///
        /// @param t Current simulation time
        /// @return Tridiagonal Hamiltonian matrix of size Dim x Dim
        ///
        /// This constructs H(t) in the qudit energy eigenbasis with:
        ///
        ///  - Diagonal energies ⟨i|H|i⟩ = E_i
        ///  - Off-diagonal dipole couplings between levels |a⟩ ↔ |b⟩
        ///  - Phase shift for Z rotations
        constexpr tridiagonal_matrix_t<Dim>
            operator()(real_t t) const noexcept
        {
            tridiagonal_matrix_t<Dim> H{};

            /* -------- Energiák -------- */
            for (natural_t i = 0; i < Dim; i++)
            {
                H[MAINDIAGONAL][i] =
                    complex_t::fromReal(m_Energies[i]);
            }

            /* -------- drive frequency -------- */
            const real_t omega_d =
                (m_Energies[m_b] - m_Energies[m_a]) + m_Detuning;

            const real_t drive =
                m_Omega * cos(omega_d * t);

            /* -------- dipole ladder -------- */
            const natural_t start = std::min(m_a, m_b);
            const natural_t end = std::max(m_a, m_b);

            for (natural_t i = start; i < end; i++)
            {
                complex_t coupling;

                switch (m_axis)
                {
                case RotationAxis::X:
                    coupling = complex_t::fromReal(drive);
                    break;

                case RotationAxis::Y:
                    coupling = complex_t(0.0, drive);
                    break;

                case RotationAxis::Z:
                    coupling = complex_t(0.0);
                    break;
                }

                H[SUPERDIAGONAL][i] = coupling;
                H[SUBDIAGONAL][i + 1] = coupling.conj();
            }

            /* -------- Z rotation -------- */
            if (m_axis == RotationAxis::Z)
            {
                H[MAINDIAGONAL][m_a] += complex_t::fromReal(-m_Omega * 0.5);
                H[MAINDIAGONAL][m_b] += complex_t::fromReal(m_Omega * 0.5);
            }

            return H;
        }
    };

}