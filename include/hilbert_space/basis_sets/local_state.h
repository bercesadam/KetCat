#pragma once
#include "core_types.h"

namespace KetCat
{
    /// @brief  Density matrix for a single local qudit subsystem.
    ///
    /// @tparam Dim  Local Hilbert-space dimension (d).
    ///
    /// @details
    /// Stores the d×d reduced density matrix ρ obtained by tracing out
    /// all qudits except the one of interest.
    template<natural_t Dim>
    struct DensityMatrix
    {
        std::array<std::array<complex_t, Dim>, Dim> m;

        constexpr void setZero() noexcept
        {
            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                    m[i][j] = complex_t::zero();
        }

        /// @brief Trace of the density matrix: Tr(ρ).
        constexpr complex_t trace() const noexcept
        {
            complex_t Tr = complex_t::zero();
            for (natural_t i = 0; i < Dim; ++i)
            {
                Tr += m[i][i];
            }
            return Tr;
        }
    };

    /// @brief  Purity of a density matrix: Tr(ρ²) ∈ (0, 1].
    ///
    /// @tparam Dim  Local Hilbert-space dimension.
    /// @param  rho  Reduced density matrix.
    /// @return      Tr(ρ²): equals 1 for a pure state, < 1 for a mixed (entangled) state.
    template<natural_t Dim>
    static constexpr real_t purity(const DensityMatrix<Dim>& rho) noexcept
    {
        real_t P = real_t(0);
        for (natural_t i = 0; i < Dim; ++i)
            for (natural_t j = 0; j < Dim; ++j)
                P += (rho.m[i][j] * rho.m[j][i]).re;   // Tr(ρ²) = Σ_ij ρ_ij ρ_ji
        return P;
    }

    /// @brief  Entanglement classification of a single qudit's reduced state.
    enum class LocalStateQualifier : natural_t
    {
        Pure      = 0,   ///< Tr(ρ²) ≈ 1  →  qudit is unentangled, has a well-defined state vector.
        Entangled = 1,   ///< Tr(ρ²) < 1  →  qudit is entangled with the rest of the register.
    };

    /// @brief  Complete local-state descriptor returned by extractLocalState.
    ///
    /// @tparam LocalDim  Local Hilbert-space dimension (d).
    ///
    /// @details
    /// Always contains the reduced density matrix ρ and its purity Tr(ρ²).
    /// When the state is pure, `pureStateVector` holds the normalized
    /// representative ket |ψ⟩ such that ρ ≈ |ψ⟩⟨ψ|.
    /// When the state is entangled, `pureStateVector` is zeroed — the reduced
    /// state cannot be described by a single ket.
    template<natural_t LocalDim>
    struct LocalStateInfo
    {
        DensityMatrix<LocalDim>            rho;             ///< Reduced density matrix (always valid).
        real_t                             purityValue;     ///< Tr(ρ²) ∈ (0, 1].
        LocalStateQualifier                     kind;            ///< Pure or Entangled classification.

        /// @brief  Representative ket for a pure local state.
        /// @note   Only meaningful when kind == LocalStateQualifier::Pure.
        StateVector<FiniteHilbertSpace<LocalDim>> pureStateVector;
    };
}
