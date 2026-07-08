#pragma once
#include <array>
#include "core_types.h"
#include "laser/laser_pulse.h"


namespace KetCat
{
    /// @brief Manages the Rotating Wave Approximation (RWA) reference frames for single and multi-atom systems.
    /// @details Calculates the static energy shifts induced by transforming the atomic Hamiltonian into the 
    /// rotating frame of the driving laser fields, tracking the relative energy offsets across multi-level manifolds.
    /// @tparam LevelCount The number of single-atom atomic energy levels.
    template<natural_t LevelCount>
    class RwaFrame
    {
    private:
        /// @brief Calculated effective stationary single-atom RWA frame energies.
        std::array<real_t, LevelCount> m_SingleRwaEnergies{};

    public:
        /// @brief Constructs the RWA frame configuration based on bare single-atom energy levels.
        /// @details Maps the physical atomic energy spectrum to a rotating frame defined by a sequence of 
        /// laser frequencies ω_i. The transformation evaluates relative energy terms as E_RWA = ΔE_bare - ∑ω_laser.
        /// @param singleAtomEnergies Array containing the bare atomic eigenenergies (e.g., Hartree fields).
        constexpr RwaFrame(const eigenenergies_t<LevelCount>& singleAtomEnergies) noexcept
        {
            real_t CumulativeOmega = 0.0;
            std::array<real_t, LevelCount - 1> LaserOmegas{};

            // Currently hardcoded for testing, TODO clean up and generalise
            LaserOmegas[0] = (singleAtomEnergies[1] - singleAtomEnergies[0]) + Units::omegaAuFromHz(1500e6);
            LaserOmegas[1] = (singleAtomEnergies[2] - singleAtomEnergies[1]) - Units::omegaAuFromHz(1500e6);
            LaserOmegas[2] = (singleAtomEnergies[3] - singleAtomEnergies[2]) + Units::omegaAuFromHz(5000e6);
            LaserOmegas[3] = (singleAtomEnergies[4] - singleAtomEnergies[3]) - Units::omegaAuFromHz(5000e6);

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                if (i > 0)
                {
                    CumulativeOmega += LaserOmegas[i - 1];
                }

                const real_t relativeEnergy = singleAtomEnergies[i] - singleAtomEnergies[0];
                m_SingleRwaEnergies[i] = relativeEnergy - CumulativeOmega;
            }
        }

        /// @brief Generates the global multi-atom RWA frame diagonal matrix representation.
        /// @details Constructs the joint uncompressed energy spectrum in the full product Hilbert space 
        /// via symmetric additive combinations of single-atom RWA frame components. 
        /// The resulting dimension scales as Dim = LevelCount^QubitCount.
        /// @tparam QubitCount Total number of physical atoms layouted in the quantum registry.
        /// @return An array populated with the joint product-space diagonal RWA energy levels.
        template<natural_t QubitCount>
        constexpr auto generateGlobalRwaEnergies() const noexcept
        {
            constexpr natural_t GlobalDim =
                ConstexprMath::pow(LevelCount, QubitCount);

            std::array<real_t, GlobalDim> GlobalRwaEnergies{};

            for (natural_t globalIdx = 0; globalIdx < GlobalDim; ++globalIdx)
            {
                real_t CurrentStateRwaEnergy = 0.0;
                natural_t tempIdx = globalIdx;

                for (natural_t q = 0; q < QubitCount; ++q)
                {
                    CurrentStateRwaEnergy += m_SingleRwaEnergies[tempIdx % LevelCount];
                    tempIdx /= LevelCount;
                }

                GlobalRwaEnergies[globalIdx] = CurrentStateRwaEnergy;
            }

            return GlobalRwaEnergies;
        }

        /// @brief Retrieve the single-atom effective RWA energy array.
        ///
        /// @return 
        ///   A reference to the calculated stationary RWA diagonals.
        constexpr const std::array<real_t, LevelCount>& getSingleRwaEnergies() const noexcept
        {
            return m_SingleRwaEnergies;
        }
    };
}