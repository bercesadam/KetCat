#pragma once
#include "core_types.h"


namespace KetCat
{
	/// @brief Master clock for synchronizing time evolution across components.
	class TimeMaster
	{
		bool m_isInitialized = false;
		real_t m_dt = 0.1;// Time step for evolution (in atomic units)

		real_t m_globalTime = 0.0;
		real_t m_currentInstructionTime = 0.0;

		// Prevent copying and moving to ensure a single authoritative time source
		TimeMaster() = default;
		~TimeMaster() = default;
		TimeMaster(const TimeMaster&) = delete;
		TimeMaster& operator=(const TimeMaster&) = delete;

	public:
		static TimeMaster& Clock()
		{
			static TimeMaster instance;
			return instance;
		}

		void init(real_t dt)
		{
			if (m_isInitialized)
			{
				return;
			}
			m_dt = dt;
			m_isInitialized = true;
		}

		real_t getTimeStep() const { return m_dt; }
		real_t getGlobalTime() const { return m_globalTime; }
		real_t getCurrentInstructionTime() const { return m_currentInstructionTime; }

		void resetCurrentInstructionClock()
		{
			m_currentInstructionTime = 0.0;
		}

		void tick()
		{
			m_globalTime += m_dt;
			m_currentInstructionTime += m_dt;
		}
	};
}
