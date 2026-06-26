#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Global simulation clock for synchronized  time evolution.
    ///
    /// @details
    ///   This singleton class provides a centralized timing source for all
    ///   time-dependent simulation components within the framework.
    ///
    ///   Responsibilities:
    ///
    ///     - Maintain the global simulation time
    ///     - Track per-instruction evolution time
    class TimeMaster
    {
        /// @brief Indicates whether the clock has been initialized.
        bool m_isInstructionStart = true;

        /// @brief Fixed TDSE evolution timestep in atomic units.
        real_t m_dt = 0.1;

        /// @brief Total accumulated simulation time.
        ///
        /// @details
        ///   Monotonically increasing global simulation clock:
        real_t m_globalTime = 0.0;

        /// @brief Elapsed time for the current logical instruction.
        ///
        /// @details
        ///   Reset after completion of each pulse command 
        real_t m_currentInstructionTime = 0.0;

        /// @brief Private constructor for singleton enforcement.
        TimeMaster() = default;

        /// @brief Private destructor for singleton enforcement.
        ~TimeMaster() = default;

        /// @brief Copy construction disabled.
        TimeMaster(const TimeMaster&) = delete;

        /// @brief Copy assignment disabled.
        TimeMaster& operator=(const TimeMaster&) = delete;

    public:
        /// @brief Access the global simulation clock instance.
        ///
        /// @return
        ///   Reference to the singleton clock instance.
        ///
        /// @details
        ///   This is the single authoritative timing source shared
        ///   by all time-dependent simulation systems.
        static TimeMaster& Clock()
        {
            static TimeMaster instance;
            return instance;
        }

        /// @brief Set the simulation timestep.
        ///
        /// @param dt
        ///   The new timestep value.
        void setTimeStep(real_t dt)
        {
            m_dt = dt;
        }

        /// @brief Retrieve the fixed simulation timestep.
        ///
        /// @return
        ///   Evolution timestep Δt in atomic units.
        real_t getTimeStep() const
        {
            return m_dt;
        }

        /// @brief Retrieve total accumulated simulation time.
        ///
        /// @return
        ///   Global simulation time.
        real_t getGlobalTime() const
        {
            return m_globalTime;
        }

        /// @brief Retrieve elapsed time for the current instruction.
        ///
        /// @return
        ///   Current instruction-local evolution time.
        real_t getCurrentInstructionTime() const
        {
            return m_currentInstructionTime;
        }

        /// @brief Determine if we are currently at the beginning of a new instruction.
        ///
        /// @return
        ///   Current instruction-local evolution time.
        bool isInstructionStart() const
        {
            return m_isInstructionStart;
        }

        /// @brief Reset the instruction-local timer.
        void resetCurrentInstructionClock()
        {
            m_currentInstructionTime = 0.0;
            m_isInstructionStart = true;
        }

        /// @brief Reset all timers; for diagnostics only.
        void reset()
        {
            m_currentInstructionTime = 0.0;
            m_globalTime = 0.0;
            m_isInstructionStart = true;
        }

        /// @brief Advance the global simulation clock by one timestep.
        void tick()
        {
            m_globalTime += m_dt;
            m_currentInstructionTime += m_dt;
            m_isInstructionStart = false;
        }
    };
}