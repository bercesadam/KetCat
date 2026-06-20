#pragma once
#include <array>
#include "core_types.h"
#include "laser/laser_pulse.h"


namespace KetCat
{
    template<natural_t LevelCount>
    class RwaFrame
    {
    private:
        std::array<real_t, LevelCount> m_SingleRwaEnergies{};

    public:
        constexpr RwaFrame(const eigenenergies_t<LevelCount>& singleAtomEnergies) noexcept
        {
            real_t CumulativeOmega = 0.0;
            std::array<real_t, LevelCount - 1> LaserOmegas{};

            // Currently hardcoded for testing, TODO clean up and generalise
            LaserOmegas[0] = (singleAtomEnergies[1] - singleAtomEnergies[0]) + Units::omegaAuFromHz(500e6);
            LaserOmegas[1] = (singleAtomEnergies[2] - singleAtomEnergies[1]) - Units::omegaAuFromHz(500e6);
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