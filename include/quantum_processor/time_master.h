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

        /// @brief Local cyclic STIRAP/Raman envelope time.
        ///
        /// @details
        ///   Used for evaluating periodic pulse envelopes:
        ///
        ///     t_cycle ∈ [0, T_cycle)
        ///
        ///   Automatically wraps once the configured instruction
        ///   time limit is reached.
        real_t m_currentStirapCycleTime = 0.0;

        /// @brief Indicates whether the current instruction timing was configured.
        bool m_isInstructionTimeInitialized = false;

        /// @brief Duration of a single pulse-cycle window.
        ///
        /// @details
        ///   Defines the wrapping limit for:
        ///     m_currentStirapCycleTime
        real_t m_oneInstructionTimeLimit = 0.0;

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

        /// @brief Set the active instruction cycle duration.
        ///
        /// @param timeLimit
        ///   Pulse or instruction duration in atomic units.
        ///
        /// @details
        ///   Defines the periodic wrapping interval used for:
        ///
        ///     m_currentStirapCycleTime
        ///
        ///   Once configured, repeated calls are ignored until
        ///   the instruction state is reset externally.
        void setInstructionTimeLimit(real_t timeLimit)
        {
            if (m_isInstructionTimeInitialized)
            {
                return;
            }

            m_oneInstructionTimeLimit = timeLimit;
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

        /// @brief Retrieve the current STIRAP cycle time.
        ///
        /// @return
        ///   Local cyclic envelope time.
        real_t getCurrentStirapCycleTime() const
        {
            return m_currentStirapCycleTime;
        }

        /// @brief Reset the instruction-local timer.
        ///
        ///   Without affecting:
        ///
        ///     • global simulation time
        ///     • STIRAP cycle timing
        void resetCurrentInstructionClock()
        {
            m_currentInstructionTime = 0.0;
        }

        /// @brief Advance the global simulation clock by one timestep.
        ///
        ///   For:
        ///
        ///     • global time
        ///     • instruction-local time
        ///     • STIRAP cycle time
        ///
        ///   STIRAP cycle wrapping:
        ///
        ///     if t_cycle ≥ T_cycle:
        ///         t_cycle → 0
        void tick()
        {
            m_globalTime += m_dt;
            m_currentInstructionTime += m_dt;
            m_currentStirapCycleTime += m_dt;

            if (m_currentStirapCycleTime >= m_oneInstructionTimeLimit)
            {
                m_currentStirapCycleTime = 0.0;
            }
        }
    };
}