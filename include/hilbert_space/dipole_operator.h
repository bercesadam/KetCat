#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"


namespace KetCat
{
    /// @brief Computes the full dipole matrix D_ij = ⟨ψ_i| x_axis |ψ_j⟩ for a basis set.
    /// 
    /// @param axis_index The spatial axis (0=x, 1=y, 2=z).
    /// @return A square matrix of size NumStates x NumStates.
    template<hilbert_space_t HilbertSpace, natural_t NumStates>
    constexpr matrix_t<NumStates> buildDipoleMatrix(const basis_set_t<HilbertSpace, NumStates>& basisStates, natural_t axis_index) noexcept
        requires (spatial_hilbert_space_t<HilbertSpace>)
    {
        constexpr natural_t GridSize = HilbertSpace::Dim;
        constexpr natural_t Steps = HilbertSpace::Steps;
        constexpr real_t dx = HilbertSpace::dx();
        const real_t dV = HilbertSpace::cellVolume(0);

        // --- Optimization: Pre-compute physical coordinates for the axis ---
        std::array<real_t, Steps> CoordinateLookup{};
        const real_t GridCenterIndex = static_cast<real_t>(Steps - 1) / 2.0;
        for (natural_t s = 0; s < Steps; ++s)
        {
            CoordinateLookup[s] = (static_cast<real_t>(s) - GridCenterIndex) * dx;
        }

        // --- Optimization: Calculate index stride ---
        natural_t Stride = 1;
        for (natural_t d = 0; d < axis_index; ++d)
        {
            Stride *= Steps;
        }

        matrix_t<NumStates> DipoleMatrix = {};

        // --- Compute matrix elements ---
        // Note: The dipole operator is Hermitian, so D_ji = conj(D_ij). 
        // We can optimize by computing only the upper triangle if needed, 
        // but for clarity, we compute the full transition matrix here.
        for (natural_t row = 0; row < NumStates; ++row)
        {
            for (natural_t col = 0; col < NumStates; ++col)
            {
                complex_t Integral = complex_t::zero();
                const auto& Bra = basisStates[row].m_Psi;
                const auto& Ket = basisStates[col].m_Psi;

                for (natural_t i = 0; i < GridSize; ++i)
                {
                    const natural_t AxisIdx = (i / Stride) % Steps;
                    const real_t Position = CoordinateLookup[AxisIdx];

                    // Standard dipole integrand: ψ*_i * x * ψ_j
                    Integral = Integral + Bra[i].conj() * Position * Ket[i];
                }

                DipoleMatrix[row][col] = Integral * dV;
            }
        }

        return DipoleMatrix;
    }
}
