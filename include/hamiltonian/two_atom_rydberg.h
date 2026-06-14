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
        /// This class is used for modeling non-local multi-qubit entangling operations, in particular
        /// Controlled-Phase (CZ) gates in neutral-atom quantum processors.
        ///
        /// @tparam LevelCount The number of modeled eigenstates per single atom!
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
            constexpr void updateMatrix(const tridiagonal_matrix_t<LevelCount>& atom0Hamiltonian,
                const tridiagonal_matrix_t<LevelCount>& atom1Hamiltonian)
            {
                create2DHamiltonian(atom0Hamiltonian, atom1Hamiltonian);

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
            constexpr void create2DHamiltonian(
                const tridiagonal_matrix_t<LevelCount>& atom0Hamiltonian,
                const tridiagonal_matrix_t<LevelCount>& atom1Hamiltonian) noexcept
            {
                m_hamiltonianMatrix = {};

                constexpr natural_t TotalDim = LevelCount * LevelCount;

                for (natural_t row = 0; row < TotalDim; ++row)
                {
                    const natural_t atom0 = row % LevelCount;
                    const natural_t atom1 = row / LevelCount;

                    // H0 ⊗ I + I ⊗ H1
                    m_hamiltonianMatrix[MAINDIAGONAL][row] =
                        atom0Hamiltonian[MAINDIAGONAL][atom0] +
                        atom1Hamiltonian[MAINDIAGONAL][atom1];

                    //
                    // atom0 transitions (±1)
                    //
                    if (atom0 < LevelCount - 1)
                    {
                        m_hamiltonianMatrix[SUPERDIAGONAL][row] =
                            atom0Hamiltonian[SUPERDIAGONAL][atom0];
                    }

                    if (atom0 > 0)
                    {
                        m_hamiltonianMatrix[SUBDIAGONAL][row] =
                            atom0Hamiltonian[SUBDIAGONAL][atom0];
                    }

                    //
                    // atom1 transitions (±LevelCount)
                    //
                    if (atom1 < LevelCount - 1)
                    {
                        m_hamiltonianMatrix[UPPER_FAR][row] =
                            atom1Hamiltonian[SUPERDIAGONAL][atom1];
                    }

                    if (atom1 > 0)
                    {
                        m_hamiltonianMatrix[LOWER_FAR][row] =
                            atom1Hamiltonian[SUBDIAGONAL][atom1];
                    }
                }
            }
        };
    }
