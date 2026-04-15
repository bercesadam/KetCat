#pragma once

#include <fstream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include "hilbert_space/state_vector.h"


namespace KetCat
{
    enum class ExportMode
    {
        RealImag,
        Probability
    };

    template <hilbert_space_t HilbertSpace>
    class StateVectorCsvExporter
    {
    public:
        using State = StateVector<HilbertSpace>;
        static constexpr natural_t Size = State::Size;

        StateVectorCsvExporter(
            const std::string& filename,
            ExportMode mode = ExportMode::RealImag)
            : m_Mode(mode)
        {
            m_File.open(filename);

            if (!m_File.is_open())
            {
                throw std::runtime_error("Failed to open CSV file: " + filename);
            }

            writeHeader();
        }

        ~StateVectorCsvExporter()
        {
            if (m_File.is_open())
            {
                m_File.close();
            }
        }

        /// @brief Append one simulation timestep to CSV
        void writeTimestep(real_t time, const State& state)
        {
            m_File << time;

            for (natural_t i = 0; i < Size; ++i)
            {
                if (m_Mode == ExportMode::RealImag)
                {
                    m_File << ',' << state[i].re;
                    m_File << ',' << state[i].im;
                }
                else
                {
                    m_File << ',' << state[i].normSquared();
                }
            }

            m_File << '\n';
        }

        /// @brief Force flush to disk
        void flush()
        {
            m_File.flush();
        }

    private:
        std::ofstream m_File;
        ExportMode m_Mode;

        void writeHeader()
        {
            m_File << "time";

            for (natural_t i = 0; i < Size; ++i)
            {
                if (m_Mode == ExportMode::RealImag)
                {
                    m_File << ",psi[" << i << "].re";
                    m_File << ",psi[" << i << "].im";
                }
                else
                {
                    m_File << ",P[" << i << "]";
                }
            }

            m_File << '\n';
        }
    };
}