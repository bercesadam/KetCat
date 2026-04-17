#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "wavefunction/wavefunction.h"


namespace KetCat
{
    /// @brief Helper function which computes the dipole matrix element ⟨ψ| x_axis |φ⟩
    ///
    /// Continuous definition:
    ///     ⟨ψ|x|φ⟩ = ∫ ψ*(r) · x · φ(r) dV
    ///
    /// Discrete approximation on a uniform grid:
    ///     Σ_i ψ_i* · x_i · φ_i · ΔV
    ///
    /// where:
    ///     x_i = physical coordinate of grid point i along selected axis,
    ///           centered such that the grid spans [-Extent/2, Extent/2].
    ///
    /// @tparam HilbertSpace Spatial Hilbert space (Uniform grid required)
    /// @param bra  ⟨ψ|
    /// @param ket  |φ⟩
    /// @param axis_index 0=x, 1=y, 2=z
    template<hilbert_space_t HilbertSpace>
    static complex_t dipoleElement(
        const StateVector<HilbertSpace>& bra,
        const StateVector<HilbertSpace>& ket,
        natural_t axis_index) noexcept
        requires (spatial_hilbert_space_t<HilbertSpace>)
    {
        constexpr natural_t Size = HilbertSpace::Dim;
        constexpr natural_t Steps = HilbertSpace::Steps;
        constexpr real_t dx = HilbertSpace::dx();
        
        // --- Pre-compute physical coordinates for the selected axis ---
        // This optimization avoids repeated floating-point divisions and subtractions 
        // inside the high-frequency loop.
        real_t CoordinateLookup[Steps];
        const real_t GridCenterIndex = static_cast<real_t>(Steps - 1) / 2.0;

        for (natural_t s = 0; s < Steps; ++s)
        {
            CoordinateLookup[s] = (static_cast<real_t>(s) - GridCenterIndex) * dx;
        }

        // --- Pre-calculate constant stride for the selected axis ---
        // Based on getIndex: Index = c[0]*1 + c[1]*Steps + c[2]*Steps^2 ...
        natural_t Stride = 1;
        for (natural_t d = 0; d < axis_index; ++d)
        {
            Stride *= Steps;
        }

        complex_t Result = complex_t::zero();

        // --- Main integration loop ---
        for (natural_t i = 0; i < Size; ++i)
        {
            // Extract the 1D grid index for the specific axis from the flat index
            const natural_t AxisGridIndex = (i / Stride) % Steps;
            const real_t Position = CoordinateLookup[AxisGridIndex];

            // Accumulate ψ*_i · x_i · φ_i
            Result = Result + bra[i].conj() * Position * ket[i];
        }

        // --- Final scaling by the uniform cell volume ---
        // For uniform grids, ΔV = dx^D is constant across all points.
        return (Result * HilbertSpace::cellVolume(0));
    }

	template<hilbert_space_t HilbertSpace, natural_t NumStates>
    class BasisSet
    {
        std::array<Wavefunction<HilbertSpace>, NumStates> m_basisStates;
        matrix_t<NumStates> m_dipoleMatrix;

        /// @brief Computes the full dipole matrix D_ij = ⟨ψ_i| x_axis |ψ_j⟩ for the basis set.
        /// 
        /// @param axis_index The spatial axis (0=x, 1=y, 2=z).
        /// @return A square matrix of size NumStates x NumStates.
        constexpr void buildDipoleMatrix(natural_t axis_index) const noexcept
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

            m_dipoleMatrix = {};

            // --- Compute matrix elements ---
            // Note: The dipole operator is Hermitian, so D_ji = conj(D_ij). 
            // We can optimize by computing only the upper triangle if needed, 
            // but for clarity, we compute the full transition matrix here.
            for (natural_t row = 0; row < NumStates; ++row)
            {
                for (natural_t col = 0; col < NumStates; ++col)
                {
                    complex_t Integral = complex_t::zero();
                    const auto& bra = m_basisStates[row].m_Psi;
                    const auto& ket = m_basisStates[col].m_Psi;

                    for (natural_t i = 0; i < GridSize; ++i)
                    {
                        const natural_t AxisIdx = (i / Stride) % Steps;
                        const real_t Position = CoordinateLookup[AxisIdx];

                        // Standard dipole integrand: ψ*_i * x * ψ_j
                        Integral = Integral + bra[i].conj() * Position * ket[i];
                    }

                    m_dipoleMatrix[row][col] = Integral * dV;
                }
            }
        }
        
    public:
        constexpr BasisSet(std::array<Wavefunction<HilbertSpace>, NumStates> basisStates)
        : m_basisStates(basisStates)
        {
            buildDipoleMatrix();
        }

        /// @brief Access a specific wavefunction in the basis
        constexpr Wavefunction<HilbertSpace>& operator[](natural_t i)
        {
            return m_basisStates[i];
        }

        /// @brief Access a specific wavefunction in the basis
        constexpr const Wavefunction<HilbertSpace>& operator[](natural_t i) const
        {
            return m_basisStates[i];
        }

        /// @brief Extracts the raw state vectors from all wavefunctions in the basis
        constexpr std::array<StateVector<HilbertSpace>, NumStates> getStateVectors() const
        {
            std::array<StateVector<HilbertSpace>, NumStates> StateVectors;
            for (natural_t i = 0; i < NumStates; ++i)
            {
                StateVectors[i] = m_basisStates[i].m_Psi;
            }
            return StateVectors;
        }

        /// @brief Extracts the eigenenergies associated with the basis states
        constexpr std::array<real_t, NumStates> getEnergies() const
        {
            std::array<real_t, NumStates> HartreeEnergies;
            for (natural_t i = 0; i < NumStates; ++i)
            {
                HartreeEnergies[i] = m_basisStates[i].m_Energy;
            }
            return HartreeEnergies;
        }
    };
}
