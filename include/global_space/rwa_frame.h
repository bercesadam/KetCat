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
        constexpr RwaFrame(const eigenenergies_t<LevelCount>& singleAtomEnergies,
                           const TwoPhotonDrive& drive) noexcept
        {
            real_t CumulativeOmega = 0.0;
            std::array<real_t, LevelCount - 1> LaserOmegas{};

            if (drive.m_groundLevelOffset < LevelCount - 1)
            {
                LaserOmegas[drive.m_groundLevelOffset] = drive.m_pump.m_omega;
            }

            if (drive.m_groundLevelOffset + 1 < LevelCount - 1)
            {
                LaserOmegas[drive.m_groundLevelOffset + 1] = drive.m_stokes.m_omega;
            }

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                if (i > 0)
                {
                    CumulativeOmega += LaserOmegas[i - 1];
                }

                const real_t relativeEnergy = singleAtomEnergies[i] - singleAtomEnergies[drive.m_groundLevelOffset];
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
    };
}