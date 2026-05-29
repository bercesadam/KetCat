#pragma once
#include "core_types.h"
#include "matrix_utils/tridiagonal_kronecker.h"
#include "operation_space/utils/matrix.h"


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
    class TwoAtomRydbergBlockage
    {
        /// @brief Global dimension of the joint 2-atom Hilbert space (Dim = LevelCount²).
        static constexpr natural_t Dim = LevelCount * LevelCount;

        /// @brief Bare energy levels of the atomic system (e.g., in Hartree or atomic units).
        std::array<real_t, LevelCount> m_Energies{};

        /// @brief Full dipole transition matrix μᵢⱼ dictating coupling strengths.
        square_matrix_t<LevelCount> m_DipoleMatrix{};

        /// @brief Cached compact pentadiagonal matrix representation of the joint 2-atom Hamiltonian.
        pentadiagonal_matrix_t<Dim> m_hamiltonianMatrix;

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
        constexpr TwoAtomRydbergBlockage(
            const real_t atomDistance,
            const natural_t rydbergLevel,
            const std::array<real_t, LevelCount>& energies,
            const square_matrix_t<LevelCount>& dipoleMatrix)
            : m_Energies(energies),
            m_DipoleMatrix(dipoleMatrix),
            m_RydbergLevelIndex(rydbergLevel)
        {
            calculateVVanDerWaals(atomDistance, rydbergLevel);
        }

        /// @brief  Assembles and returns the joint pentadiagonal Hamiltonian matrix.
        /// @param  singleAtomRydbergExcitation Single-atom driving matrix (Rabi frequencies and detunings).
        /// @return Const reference or value of the populated pentadiagonal_matrix_t.
        ///
        /// @details
        /// Maps individual atomic operators to the tensor space via H = H₁ ⊗ I + I ⊗ H₂.
        /// The Rydberg interaction is then dynamically added as a local energy shift to the 
        /// diagonal element representing the doubly excited state:
        ///
        ///   VrrIndex = r · LevelCount + r
        const pentadiagonal_matrix_t<Dim> getMatrix(const tridiagonal_matrix_t<LevelCount>& singleAtomRydbergExcitation)
        {
            tridiagonal_matrix_t<LevelCount> I{};
            std::fill(I[MAINDIAGONAL].begin(),
                I[MAINDIAGONAL].end(),
                complex_t::fromReal(1.0));

            pentadiagonal_matrix_t<Dim> Atom1Excitation = tensorProduct(singleAtomRydbergExcitation, I);
            pentadiagonal_matrix_t<Dim> Atom2Excitation = tensorProduct(I, singleAtomRydbergExcitation);

            // Directly compute linear combination using continuous memory blocks
            for (natural_t d = 0; d < 5U; ++d)
            {
                for (natural_t i = 0; i < Dim; ++i)
                {
                    m_hamiltonianMatrix[d][i] = Atom1Excitation[d][i] + Atom2Excitation[d][i];
                }
            }

            // Apply the non-local Rydberg interaction penalty V_vdW to the diagonal state |r,r⟩
            const natural_t VrrIndex = m_RydbergLevelIndex * LevelCount + m_RydbergLevelIndex;
            m_hamiltonianMatrix[MAINDIAGONAL][VrrIndex] =
                m_hamiltonianMatrix[MAINDIAGONAL][VrrIndex] + complex_t::fromReal(m_VVanDerWaals);

            return m_hamiltonianMatrix;
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
        void calculateVVanDerWaals(real_t atomDistanceR, natural_t rydbergIndex)
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
                    if (std::abs(Denominator) < 1e-12) continue;

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
    };
}
