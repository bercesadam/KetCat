#pragma once

#include <bitset>
#include <cmath>
#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "local_space/neutral_atom_config.h"


namespace KetCat
{
    /// @brief  Extracts the computational state vector coefficients (amplitudes and phases).
    ///         Optionally renormalizes the results if probability leaked outside the computational subspace.
    /// @tparam QubitCount  Total number of qubits in the system.
    /// @tparam Config      Compile-time physical configuration of the neutral atom.
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class LogicalStateVectorExtractor
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;
        using FullStateVectorType =
            StateVector<typename GlobalStateManager::FullHilbertSpace, QuantumPicture::Schrodinger>;

        static constexpr natural_t LogicalDimSize = ConstexprMath::pow(natural_t(2), QubitCount);

    public:
        /// @brief Extracts the complex amplitudes of the logical state vector from the full global state vector.
        ///        This preserves both amplitude and phase information for phase disk visualizations.
        static constexpr std::array<complex_t, LogicalDimSize>
            extractLogicalStateVector(const FullStateVectorType& psi, bool renormalize) noexcept
        {
            std::array<complex_t, LogicalDimSize> LogicalAmplitudes{};
            real_t TotalLogicalProbability = 0.0;

            // Little-Endian mapping: logicalIdx binary representation matches qubit indexing directly
            for (natural_t logicalIdx = 0; logicalIdx < LogicalDimSize; ++logicalIdx)
            {
                const std::bitset<QubitCount> Bitstring(logicalIdx);

                natural_t GlobalIndex = 0;
                natural_t Multiplier = 1;

                for (natural_t i = 0; i < QubitCount; ++i)
                {
                    const bool Bit = Bitstring[i];
                    const natural_t PhysicalLevel = (Bit ? ConfigType::Logical1Level : ConfigType::Logical0Level);

                    GlobalIndex += PhysicalLevel * Multiplier;
                    Multiplier *= ConfigType::LevelCount;
                }

                // Extract the full complex amplitude containing the phase information
                const complex_t Amp = psi[GlobalIndex];
                LogicalAmplitudes[logicalIdx] = Amp;

                // Track total probability within the logical subspace for leakage normalization
                TotalLogicalProbability += Amp.normSquared();
            }

            // Phase-preserving renormalization to correct for physical leakage errors
            if (renormalize && TotalLogicalProbability > 0.0)
            {
                // Since total population inside the computational subspace has decreased,
                // amplitudes must be scaled by the square root of the total remaining probability.
                // This correctly scales the magnitudes while leaving the complex phases intact.
                const real_t NormFactor = 1.0 / std::sqrt(TotalLogicalProbability);
                for (natural_t i = 0; i < LogicalDimSize; ++i)
                {
                    LogicalAmplitudes[i] = LogicalAmplitudes[i] * NormFactor;
                }
            }

            return LogicalAmplitudes;
        }
    };
}
