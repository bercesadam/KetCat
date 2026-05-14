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
        bool m_isInitialized = false;

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

        /// @brief Initialize the simulation timestep.
        ///
        /// @param dt
        ///   Fixed evolution timestep in atomic units.
        ///
        /// @details
        ///   Initialization is performed only once.
        ///
        ///   Subsequent calls are ignored in order to preserve
        ///   global timing consistency across the simulation.
        ///
        void init(real_t dt)
        {
            if (m_isInitialized)
            {
                return;
            }

            m_dt = dt;
            m_isInitialized = true;
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

        /// @brief Reset the instruction-local timer.
        void resetCurrentInstructionClock()
        {
            m_currentInstructionTime = 0.0;
        }

        /// @brief Advance the global simulation clock by one timestep.
        void tick()
        {
            m_globalTime += m_dt;
            m_currentInstructionTime += m_dt;
        }
    };
}