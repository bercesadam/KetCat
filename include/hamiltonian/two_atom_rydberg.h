#pragma once 
#include "core_types.h"
#include "matrix_utils/tridiagonal_kronecker.h"
#include "operation_space/utils/matrix.h"


namespace KetCat
{
    

    template <natural_t LevelCount>
    class TwoAtomRydbergBlockage
    {
        static constexpr natural_t Dim = LevelCount * LevelCount;

        /// @brief Bare energy levels of the system (e.g. Hartree units).
        std::array<real_t, LevelCount> m_Energies{};

        /// @brief Full dipole transition matrix μᵢⱼ.
        square_matrix_t<LevelCount> m_DipoleMatrix{};

        /// @brief Two-atom Hamiltonian matrix
        Matrix<Dim> m_hamiltonianMatrix;

        ///@brief
        natural_t m_RydbergLevelIndex;

        /// @brief 
        real_t m_VVanDerWaals;

    public:
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
            
        const Matrix<Dim> getMatrix(const tridiagonal_matrix_t<LevelCount>& singleAtomRydbergExcitation)
        {
            static const tridiagonal_matrix_t<LevelCount> I = []() {
                tridiagonal_matrix_t<LevelCount> tmp{};
                std::fill(tmp[MAINDIAGONAL].begin(),
                    tmp[MAINDIAGONAL].end(),
                    complex_t::fromReal(1.0));
                return tmp;
            }();

            Matrix<Dim> Atom1Excitation = tensorProduct(singleAtomRydbergExcitation, I);
            Matrix<Dim> Atom2Excitation = tensorProduct(I, singleAtomRydbergExcitation);

            constexpr natural_t TotalElements = Dim * Dim;
            const auto* Ptr1 = &Atom1Excitation.m[0][0];
            const auto* Ptr2 = &Atom2Excitation.m[0][0];
            auto* Dest = &m_hamiltonianMatrix.m[0][0];

            for (natural_t i = 0; i < TotalElements; ++i)
            {
                Dest[i] = Ptr1[i] + Ptr2[i];
            }

			static const natural_t VrrIndex = m_RydbergLevelIndex * LevelCount + m_RydbergLevelIndex;
			m_hamiltonianMatrix.at(VrrIndex, VrrIndex) =
                m_hamiltonianMatrix.at(VrrIndex, VrrIndex) + complex_t::fromReal(m_VVanDerWaals);

			return m_hamiltonianMatrix;
        }

    private:
        void calculateVVanDerWaals(real_t atomDistanceR, natural_t rydbergIndex)
        {
            real_t Sum = 0.0;
            real_t E_r = m_Energies[rydbergIndex];

            // Végigmegyünk az összes lehetséges kétatomos (n1, n2) kombináción
            for (natural_t n1 = 0; n1 < LevelCount; ++n1)
            {
                for (natural_t n2 = 0; n2 < LevelCount; ++n2)
                {
                    // Kihagyjuk a kiindulási állapotot, hogy ne legyen nullával osztás
                    if (n1 == rydbergIndex && n2 == rydbergIndex) continue;

                    real_t E_n1 = m_Energies[n1];
                    real_t E_n2 = m_Energies[n2];
                    real_t Denominator = (2.0 * E_r) - (E_n1 + E_n2);

                    // Ha az energiakülönbség elhanyagolható, a perturbációszámítás szingulárissá válna
                    if (std::abs(Denominator) < 1e-12) continue;

                    // Dipól mátrixelemek lekérése (r -> n1 és r -> n2)
                    real_t mu_1 = m_DipoleMatrix[rydbergIndex][n1].re;
                    real_t mu_2 = m_DipoleMatrix[rydbergIndex][n2].re;

                    // Geometriai szorzó (ha z-tengely menti a csatolás, a faktor 4.0)
                    real_t matrixElementSq = 4.0 * (mu_1 * mu_1) * (mu_2 * mu_2);

                    Sum += matrixElementSq / Denominator;
                }
            }

            // V_vdw = C6 / R^6
            real_t r6 = std::pow(atomDistanceR, 6.0);
            m_VVanDerWaals = Sum / r6;
        }
    };
}