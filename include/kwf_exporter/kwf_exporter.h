#pragma once
#include <fstream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <iostream>

#include "simulation_data.h"
#include "hilbert_space/state_vector.h"



/// @file
/// This approach replaces the former outputs which were based on visualizations
/// generated directly by KetCat executables. All outputs can be visualized by external 
/// Python-based tools which are more suitable and faster for visualization purposes.
///
/// KetCat Binary Wavefunction Format  (.kwf)
/// ──────────────────────────────────────────
/// Strings (caption) are stored as plain ASCII/UTF-8 with a 2-byte
/// length prefix.  All numeric data (time, amplitudes, probabilities)
/// are stored as raw IEEE-754 little-endian binary, eliminating the
/// decimal-conversion overhead of CSV and cutting file size by ~3-4x.
///
/// Layout
/// ──────
///   [FILE HEADER]
///     4 bytes  – magic bytes "KWF\x02"  (format version 2)
///     1 byte   – ExportMode: 0 = RealImag, 1 = Probability
///     1 byte   – uint8  : number of qubits in the experiment
///     8 bytes  – uint64 : number of basis states (Size)
///     8 bytes  – uint64 : Spatial discretization steps of the underlying Hilbert space
///     8 bytes  – double : Physical extent of the underlying Hilbert space in AU
///
///   [FRAME]  (repeated once per writeTimestep call)
///     2 bytes  – uint16 : byte-length N of the caption string
///     N bytes  – UTF-8  : caption (no null terminator)
///     8 bytes  – float64: simulation time
///     Per quantum state state_dim: 2 x float64 (re, im) [Global State Vector]
///     Per qubit num_qubits:        1 x float64 (purity)
///     Laser params:                4 x float64 (w1, i1, w2, i2)
///     Per basis state, RealImag mode:  2 x float64  (re, im)
///     Per basis state, Probability mode: 1 x float64 (|psi|^2)
///
/// All multi-byte integers and floats are little-endian.

namespace KetCat
{
    /// @brief Selects what wavefunction data is written per frame.
    enum class ExportMode
    {
        RealImag,
        Probability
    };

    /// @brief Serialises StateVector snapshots to a compact binary file (.kwf).
    /// @tparam HilbertSpace  The Hilbert-space type tag; must satisfy hilbert_space_t.
    template <hilbert_space_t HilbertSpace, natural_t QubitCount>
    class StateVectorExporter
    {
		/// @brief Binary file stream for writing the .kwf output.
        std::ofstream m_File;

		/// @brief Export mode determining the per-state payload format (Real+Imaginary or Probability).
        ExportMode    m_Mode;

		/// @brief Number of qubits in the simulation, needed to determine the per-frame payload size.
        uint8_t       m_NumQubits;

    public:
        using State = StateVector<HilbertSpace>;

        /// @brief Total number of basis states in the Hilbert space.
        static constexpr natural_t Size = State::Size;

        /// @brief Opens (or creates) a .kwf file and writes the file header.
        StateVectorExporter(
            const std::string& filename,
            uint8_t numQubits,
            ExportMode mode = ExportMode::RealImag)
            : m_Mode(mode), m_NumQubits(numQubits)
        {
            m_File.open(filename, std::ios::binary | std::ios::out);

            if (!m_File.is_open())
            {
                std::cout << "Failed to open binary export file: " + filename << std::endl;
                exit(1);
            }

            writeFileHeader();
        }

        /// @brief Closes the file on destruction.
        ~StateVectorExporter()
        {
            if (m_File.is_open())
                m_File.close();
        }

        /// @brief Appends one simulation timestep to the binary file.
        void writeTimestep(const SimulationView<HilbertSpace, QubitCount>& view)
        {
            // Simulation time as a 64-bit IEEE-754 double
            const double t = static_cast<double>(view.m_time);
            writeRaw(t);

            // Write full global state vector amplitudes
            const size_t stateDim = view.m_outputProbabilities.size();
            for (size_t i = 0; i < stateDim; ++i)
            {
                const double p = static_cast<double>(view.m_outputProbabilities[i]);
                writeRaw(p);
            }

			// Write captions for each qubit
            for (size_t i = 0; i < m_NumQubits; ++i)
            {
                // Caption: 2-byte length prefix + raw UTF-8 bytes (no null terminator)
                const auto captionLen = static_cast<uint16_t>(view.m_qubitDatum[i].m_title.size());
                writeRaw(captionLen);
                m_File.write(view.m_qubitDatum[i].m_title.data(), captionLen);
            }

            // Write Bloch vectors for each qubit
            for (size_t i = 0; i < m_NumQubits; ++i)
            {
                writeRaw(static_cast<double>(view.m_qubitDatum[i].m_alpha.re));
                writeRaw(static_cast<double>(view.m_qubitDatum[i].m_alpha.im));
                writeRaw(static_cast<double>(view.m_qubitDatum[i].m_beta.re));
                writeRaw(static_cast<double>(view.m_qubitDatum[i].m_beta.im));
            }

            // Write purity values for each qubit
            for (size_t i = 0; i < m_NumQubits; ++i)
            {
                const double purity = static_cast<double>(view.m_qubitDatum[i].m_purity);
                writeRaw(purity);
            }

            // Write laser parameters for each individual atom
            for (size_t i = 0; i < m_NumQubits; ++i)
            {
                const double l1w = static_cast<double>(view.m_qubitDatum[i].m_laser1Wavelength);
                const double l1i = static_cast<double>(view.m_qubitDatum[i].m_laser1Intensity);
                const double l2w = static_cast<double>(view.m_qubitDatum[i].m_laser2Wavelength);
                const double l2i = static_cast<double>(view.m_qubitDatum[i].m_laser2Intensity);
                writeRaw(l1w);
                writeRaw(l1i);
                writeRaw(l2w);
                writeRaw(l2i);
            }

            // Wavefunction payload written sequentially for each atom
            for (uint8_t q = 0; q < m_NumQubits; ++q)
            {
                for (natural_t i = 0; i < Size; ++i)
                {
                    if (m_Mode == ExportMode::RealImag)
                    {
                        const double re = static_cast<double>(view.m_qubitDatum[q].m_psi2D[i].re);
                        const double im = static_cast<double>(view.m_qubitDatum[q].m_psi2D[i].im);
                        writeRaw(re);
                        writeRaw(im);
                    }
                    else
                    {
                        const double p = static_cast<double>(view.m_qubitDatum[q].m_psi2D[i].normSquared());
                        writeRaw(p);
                    }
                }
            }
        }

        /// @brief Forces all buffered frame data to be written to disk.
        void flush()
        {
            m_File.flush();
        }

    private:
        /// @brief Writes any trivially-copyable value as raw little-endian bytes.
        template <typename T>
        void writeRaw(const T& value)
        {
            static_assert(std::is_trivially_copyable_v<T>);
            m_File.write(reinterpret_cast<const char*>(&value), sizeof(T));
        }

        /// @brief Writes the file header exactly once at construction time.
        void writeFileHeader()
        {
            // Magic (v2.0)
            m_File.write("KWF\x02", 4);

            // Single-byte mode flag so the reader knows the per-state payload size
            const uint8_t modeByte = (m_Mode == ExportMode::RealImag) ? 0u : 1u;
            writeRaw(modeByte);

            // Number of qubits in the simulation
            writeRaw(m_NumQubits);

            // Basis-state count lets the reader size its buffers
            const uint64_t size = static_cast<uint64_t>(Size);
            writeRaw(size);

            const uint64_t steps = static_cast<uint64_t>(HilbertSpace::Steps);
            writeRaw(steps);

            const double extent = static_cast<double>(HilbertSpace::Extent);
            writeRaw(extent);
        }
    };
}