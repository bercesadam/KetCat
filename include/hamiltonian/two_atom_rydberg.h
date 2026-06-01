#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Hamiltonian builder for a interacting two-atom system with Rydberg blockage.
    ///
    /// @details
    /// This class models the coherent interaction of two multi-level atoms under the influence 
    /// of a long-range van der Waals interaction. The primary physical effect is the Rydberg blockade, 
    /// where the simultaneous excitation of both atoms to a high-lying Rydberg state |r,r⟩ is shifted 
    /// out of resonance by an energy penalty V_vdW.
    ///
    /// The joint two-atom Hamiltonian is constructed in the product basis |i,j⟩ = |i⟩ ⊗ |j⟩ using 
    /// Kronecker tensor products of single-atom operators. Thanks to the decoupled nature of the 
    /// independent 1D atomic excitations, the non-zero elements are rigorously confined to a sávos 
    /// structure, allowing the use of the highly optimized pentadiagonal framework:
    ///
    ///   H_tot = (H_atom1 ⊗ I) + (I ⊗ H_atom2) + V_vdW · |r,r⟩⟨r,r|
    ///
    /// This solver is crucial for modeling non-local multi-qubit entangling operations, such as 
    /// Controlled-Phase (CZ) gates in neutral-atom quantum processors.
    ///
    /// @tparam LevelCount The number of internal electronic energy levels resolved per single atom.
    template <natural_t LevelCount>
    class TwoAtomRydbergBlockade
    {
        /// @brief Global dimension of the joint 2-atom Hilbert space (LevelCount = LevelCount²).
        static constexpr natural_t OutputDim = LevelCount * LevelCount;

        /// @brief Bare energy levels of the atomic system (e.g., in Hartree or atomic units).
        std::array<real_t, LevelCount> m_Energies{};

        /// @brief Full dipole transition matrix μᵢⱼ dictating coupling strengths.
        square_matrix_t<LevelCount> m_DipoleMatrix{};

        /// @brief Cached compact pentadiagonal matrix representation of the joint 2-atom Hamiltonian.
        five_band_matrix_t<LevelCount> m_hamiltonianMatrix;

        /// @brief Target electronic index corresponding to the highly excited Rydberg state |r⟩.
        natural_t m_RydbergLevelIndex;

        /// @brief Computed van der Waals interaction energy shift V_vdW = C₆ / R⁶.
        real_t m_VVanDerWaals;

    public:
        /// @brief  Constructs a two-atom Rydberg system configuration.
        /// @param  atomDistance   Interatomic spatial distance R.
        /// @param  rydbergLevel   The state index designating the target Rydberg level.
        /// @param  energies       The bare single-atom energy spectrum vector.
        /// @param  dipoleMatrix   The transition dipole moment matrix.
        constexpr TwoAtomRydbergBlockade(
            const real_t atomDistance,
            const natural_t rydbergLevel,
            const std::array<real_t, LevelCount>& energies,
            const square_matrix_t<LevelCount>& dipoleMatrix)
            : m_Energies(energies),
            m_DipoleMatrix(dipoleMatrix),
            m_RydbergLevelIndex(rydbergLevel)
        {
            calculateVVanDerWaals(atomDistance, rydbergLevel);
			std::cout << "Rydberg blockade energy shift V_vdW: " << m_VVanDerWaals << " atomic units" << std::endl;
        }

        /// @brief Get the current Hamiltonian matrix.
        ///
        /// @return
        ///   Reference to the internal tridiagonal matrix.
        constexpr const five_band_matrix_t<LevelCount>& getMatrix() const noexcept
        {
            return m_hamiltonianMatrix;
        }

        /// @brief  Assembles and returns the joint pentadiagonal Hamiltonian matrix.
        /// @param  singleAtomRydbergExcitation Single-atom driving matrix (Rabi frequencies and detunings).
        /// @return Const reference or value of the populated five_band_matrix_t.
        ///
        /// @details
        /// Maps individual atomic operators to the tensor space via H = H₁ ⊗ I + I ⊗ H₂.
        /// The Rydberg interaction is then dynamically added as a local energy shift to the 
        /// diagonal element representing the doubly excited state:
        ///
        ///   VrrIndex = r · LevelCount + r
        constexpr void updateMatrix(const tridiagonal_matrix_t<LevelCount>& singleAtomRydbergExcitation)
        {
            create2DHamiltonian(singleAtomRydbergExcitation);

            // Apply the non-local Rydberg interaction penalty V_vdW to the diagonal state |r,r⟩
            const natural_t VrrIndex = m_RydbergLevelIndex * LevelCount + m_RydbergLevelIndex;
            m_hamiltonianMatrix[MAINDIAGONAL][VrrIndex] =
                m_hamiltonianMatrix[MAINDIAGONAL][VrrIndex] + complex_t::fromReal(m_VVanDerWaals);
        }

    private:
        /// @brief  Evaluates the van der Waals interaction shift via second-order perturbation theory.
        /// @param  atomDistanceR Interatomic separation distance R.
        /// @param  rydbergIndex  The state vector index corresponding to the target Rydberg channel.
        ///
        /// @details
        /// Computes the effective asymptotic dispersion coefficient C₆ by summing over all 
        /// virtual intermediate channels |n₁,n₂⟩ generated by dipole-dipole transitions:
        ///
        ///   C₆ = ∑  |⟨r,r| Û_dd |n₁,n₂⟩|² / (2E_r - E_n1 - E_n2)
        ///
        /// This assumes a geometric coupling coefficient of 4.0 matching quantizations along the 
        /// interatomic z-axis. The final energy penalty isscaled as V_vdW = C₆ / R⁶.
        constexpr void calculateVVanDerWaals(real_t atomDistanceR, natural_t rydbergIndex)
        {
            real_t Sum = 0.0;
            real_t E_r = m_Energies[rydbergIndex];

            // Loop over all possible virtual two-atom (n1, n2) combinations
            for (natural_t n1 = 0; n1 < LevelCount; ++n1)
            {
                for (natural_t n2 = 0; n2 < LevelCount; ++n2)
                {
                    // Skip the initial state to prevent division by zero singularities
                    if (n1 == rydbergIndex && n2 == rydbergIndex) continue;

                    real_t E_n1 = m_Energies[n1];
                    real_t E_n2 = m_Energies[n2];
                    real_t Denominator = (2.0 * E_r) - (E_n1 + E_n2);

                    // Protect against near-resonant degeneracies where perturbation theory fails
                    if (ConstexprMath::abs(Denominator) < 1e-12) continue;

                    // Fetch single-atom transition dipole matrix elements (r -> n1 and r -> n2)
                    real_t mu_1 = m_DipoleMatrix[rydbergIndex][n1].re;
                    real_t mu_2 = m_DipoleMatrix[rydbergIndex][n2].re;

                    // Angular geometric prefix factor (equals 4.0 for pure axial interactions)
                    real_t matrixElementSq = 4.0 * (mu_1 * mu_1) * (mu_2 * mu_2);

                    Sum += matrixElementSq / Denominator;
                }
            }

            // Apply scaling: V_vdw = C6 / R^6
            real_t r6 = std::pow(atomDistanceR, 6.0);
            m_VVanDerWaals = Sum / r6;
        }

        /// @brief Constructs a 2D Hamiltonian pentadiagonal matrix from a 1D tridiagonal operator
        ///        by computing the symmetric operator sum: m_hamiltonianMatrix = H_1D ⊗ I + I ⊗ H_1D
        ///
        /// @details
        /// This function exploits the exact sparse layout of the 2D Hamiltonian for two identical
        /// coupled subsystems (atoms). It maps the 1D physical couplings directly into the 
        /// optimized 5-band matrix framework in O(N²) time complexity, completely avoiding 
        /// the cross-diagonal terms that an unconstrained tensor product would generate.
        ///
        /// @tparam LevelCount Dimension of the 1D Hilbert space (number of energy levels for one atom).
        /// @param H_1D        The 1D tridiagonal Hamiltonian operator matrix.
        /// @return            The populated, compact five_band_matrix_t representation for the 2D system.
        constexpr void create2DHamiltonian(const tridiagonal_matrix_t<LevelCount>& H_1D) noexcept
        {
            m_hamiltonianMatrix = {};

            constexpr natural_t TotalDim = LevelCount * LevelCount;

            // 2. A 2D Kronecker-szorzatok felépítése O(N) sávos bejárással
            // 2. Sávok feltöltése a solver belső indexelési sémájának megfelelően
            for (natural_t Row = 0; Row < TotalDim; ++Row)
            {
                const natural_t i = Row / LevelCount; // Első atom állapota
                const natural_t n = Row % LevelCount; // Második atom állapota

                // FŐÁTLÓ: H_1D ⊗ I + I ⊗ H_1D
                m_hamiltonianMatrix[MAINDIAGONAL][Row] = H_1D[MAINDIAGONAL][i] + H_1D[MAINDIAGONAL][n];

                // KÖZELI SÁVOK: I ⊗ H_1D (Második atom ugrásai a blokkon belül)
                if (n < LevelCount - 1)
                {
                    m_hamiltonianMatrix[SUPERDIAGONAL][Row] = H_1D[SUPERDIAGONAL][n];
                }
                if (n > 0)
                {
                    // A tridiagonális mátrixod konvenciója alapján a SUBDIAGONAL[Row] 
                    // a (Row, Row-1) elemet jelenti. Mivel a blokkon belül vagyunk, n indexel.
                    m_hamiltonianMatrix[SUBDIAGONAL][Row] = H_1D[SUBDIAGONAL][n];
                }

                // TÁVOLI SÁVOK: H_1D ⊗ I (Első atom ugrásai a blokkok között)
                // Ha i < LevelCount - 1, akkor van ugrás előre (i -> i+1) blokkszinten.
                // Ez a solvernél az UPPER_FAR[Row] sávban lakik.
                if (i < LevelCount - 1)
                {
                    m_hamiltonianMatrix[UPPER_FAR][Row] = H_1D[SUPERDIAGONAL][i];
                }
                // Ha i > 0, akkor van ugrás hátra (i -> i-1) blokkszinten.
                // Ez a solvernél a LOWER_FAR[Row] sávban lakik.
                if (i > 0)
                {
                    m_hamiltonianMatrix[LOWER_FAR][Row] = H_1D[SUBDIAGONAL][i];
                }
            }
        }
    };
}
