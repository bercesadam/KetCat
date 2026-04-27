#pragma once

#include <fstream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include "hilbert_space/state_vector.h"


/// @file
/// This approach replaces the former outputs which were based on visualizations
/// generated directly by KetCat executables. ALl outputs can be visualized by external 
/// Python-based tools which are more suitable and faster for visualization purposes.
///
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
///     4 bytes  – magic bytes "KWF\x01"  (format version 1)
///     1 byte   – ExportMode: 0 = RealImag, 1 = Probability
///     8 bytes  – uint64 : number of basis states (Size)
///     8 bytes  – uint64 : Spatial discretization steps of the underlying Hilbert space
///     8 bytes  – double : Physical extent of the underlying Hilbert space in AU
///
///   [FRAME]  (repeated once per writeTimestep call)
///     2 bytes  – uint16 : byte-length N of the caption string
///     N bytes  – UTF-8  : caption (no null terminator)
///     8 bytes  – float64: simulation time
///     Per basis state, RealImag mode:  2 x float64  (re, im)
///     Per basis state, Probability mode: 1 x float64 (|psi|^2)
///
/// All multi-byte integers and floats are little-endian.
////

namespace KetCat
{
    /// @brief Selects what wavefunction data is written per frame.
    /// @details RealImag stores the full complex amplitude (re, im) for each
    ///          basis state, preserving all phase information.
    ///          Probability stores only |psi_i|^2, halving the data volume
    ///          at the cost of discarding phase.
    enum class ExportMode
    {
        RealImag,
        Probability
    };

    /// @brief Serialises StateVector snapshots to a compact binary file (.kwf).
    /// @tparam HilbertSpace  The Hilbert-space type tag; must satisfy hilbert_space_t.
    ///
    /// @details
    ///   Drop-in replacement for StateVectorCsvExporter with an identical
    ///   public interface.  Instead of human-readable CSV text, it writes a
    ///   tightly-packed binary stream where captions are stored as length-
    ///   prefixed ASCII/UTF-8 and all numeric values are raw IEEE-754
    ///   little-endian doubles.
    ///
    ///   For a 256x256 grid in RealImag mode, every frame occupies
    ///   exactly  256*256/// 2/// 8 = 1,048,576 bytes  versus the ~3-4 MB
    ///   that CSV text typically requires for the same data.
    ///
    ///   Usage:
    ///   @code
    ///   StateVectorExporter<MySpace> exporter("run.kwf");
    ///   exporter.writeTimestep(t, psi, "step 0");
    ///   exporter.flush();
    ///   @endcode
    template <hilbert_space_t HilbertSpace>
    class StateVectorExporter
    {
    public:
        using State = StateVector<HilbertSpace>;

        /// @brief Total number of basis states in the Hilbert space.
        static constexpr natural_t Size = State::Size;

        /// @brief Opens (or creates) a .kwf file and writes the file header.
        /// @param filename  Path to the output file.  The file is truncated if it
        ///                  already exists.
        /// @param mode      ExportMode::RealImag (default) or ExportMode::Probability.
        /// @throws std::runtime_error if the file cannot be opened.
        StateVectorExporter(
            const std::string& filename,
            ExportMode mode = ExportMode::RealImag)
            : m_Mode(mode)
        {
            m_File.open(filename, std::ios::binary | std::ios::out);

            if (!m_File.is_open())
                throw std::runtime_error("Failed to open binary export file: " + filename);

            writeFileHeader();
        }

        /// @brief Closes the file on destruction.
        /// @details Any buffered data is flushed automatically by the ofstream
        ///          destructor; call flush() explicitly if you need an earlier
        ///          guarantee (e.g. between timesteps in a long run).
        ~StateVectorExporter()
        {
            if (m_File.is_open())
                m_File.close();
        }

        /// @brief Appends one simulation timestep to the binary file.
        /// @param time   Simulation time corresponding to this snapshot.
        /// @param state  The StateVector to be exported.
        /// @param title  Human-readable caption stored as length-prefixed UTF-8.
        ///               Maximum caption length is 65535 bytes (uint16 limit).
        ///
        /// @details Frame binary layout:
        ///   uint16  caption_length
        ///   char[]  caption  (caption_length bytes, no null terminator)
        ///   float64 time
        ///   float64 re[0], float64 im[0], ...  (RealImag mode)
        ///   -- or --
        ///   float64 |psi[0]|^2, ...             (Probability mode)
        void writeTimestep(real_t time, const State& state, const std::string& title,
            const complex_t alpha, const complex_t beta)
        {
            // Caption: 2-byte length prefix + raw UTF-8 bytes (no null terminator)
            const auto captionLen = static_cast<uint16_t>(title.size());
            writeRaw(captionLen);
            m_File.write(title.data(), captionLen);

            // Simulation time as a 64-bit IEEE-754 double
            const double t = static_cast<double>(time);
            writeRaw(t);

            // Write Bloch vectors for Bloch sphere visu
            writeRaw(alpha.re);
			writeRaw(alpha.im);
			writeRaw(beta.re);
			writeRaw(beta.im);

            // Wavefunction payload
            for (natural_t i = 0; i < Size; ++i)
            {
                if (m_Mode == ExportMode::RealImag)
                {
                    const double re = static_cast<double>(state[i].re);
                    const double im = static_cast<double>(state[i].im);
                    writeRaw(re);
                    writeRaw(im);
                }
                else
                {
                    const double p = static_cast<double>(state[i].normSquared());
                    writeRaw(p);
                }
            }
        }

        /// @brief Forces all buffered frame data to be written to disk.
        /// @details Useful after writing a checkpoint frame during a long
        ///          simulation so the file remains valid if the process is
        ///          interrupted.
        void flush()
        {
            m_File.flush();
        }

    private:
        std::ofstream m_File;
        ExportMode    m_Mode;

        /// @brief Writes any trivially-copyable value as raw little-endian bytes.
        /// @tparam T  Must be trivially copyable (enforced by static_assert).
        /// @note On virtually all modern x86/ARM platforms, the native float
        ///       representation is already little-endian.  Add a byte-swap
        ///       step here if big-endian portability is ever required.
        template <typename T>
        void writeRaw(const T& value)
        {
            static_assert(std::is_trivially_copyable_v<T>);
            m_File.write(reinterpret_cast<const char*>(&value), sizeof(T));
        }

        /// @brief Writes the 13-byte file header exactly once at construction time.
        /// @details Header layout:
        ///   [0..3]   magic "KWF\x01"
        ///   [4]      uint8  ExportMode  (0 = RealImag, 1 = Probability)
        ///   [5..12]  uint64 Size        (number of basis states)
        ///   [13..20] uint64 Steps       (spatial discretization steps of the underlying Hilbert space)
        ///   [21..28] double Extent      (physical extent of the Hilbert space in AU)
        void writeFileHeader()
        {
            // Magic bytes identify the format and encode the format version
            m_File.write("KWF\x01", 4);

            // Single-byte mode flag so the reader knows the per-state payload size
            const uint8_t modeByte = (m_Mode == ExportMode::RealImag) ? 0u : 1u;
            writeRaw(modeByte);

            // Basis-state count lets the reader size its buffers without
            // scanning the entire file first
            const uint64_t size = static_cast<uint64_t>(Size);
            writeRaw(size);

            const uint64_t steps = static_cast<uint64_t>(HilbertSpace::Steps);
            writeRaw(steps);

            const double extent = static_cast<double>(HilbertSpace::Extent);
            writeRaw(extent);
        }
    };
}
