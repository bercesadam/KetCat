#pragma once
#include <bitset>
#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "local_space/neutral_atom_config.h"


namespace KetCat
{
    /// @brief  Extracts the measurement probabilities of the logical state vector.
    ///         Optionally renormalizes the results if probability leaked outside the computational subspace.
    /// @tparam QubitCount  Total number of qubits in the system.
    /// @tparam Config      Compile-time physical configuration of the neutral atom.
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class LogicalProbabilityExtractor
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;
        using FullStateVectorType = StateVector<typename GlobalStateManager::FullHilbertSpace>;

        static constexpr natural_t LogicalDimSize = ConstexprMath::pow(natural_t(2), QubitCount);

    public:
        /// @brief Extracts the probabilities of the logical state vector from the full global state vector.
        static constexpr probability_vector_t<LogicalDimSize>
            extractLogicalProbabilities(const FullStateVectorType& psi, bool renormalize) noexcept
        {
            probability_vector_t<LogicalDimSize> LogicalProbabilities{};
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

                const real_t P = psi[GlobalIndex].normSquared();
                LogicalProbabilities[logicalIdx] = P;
                TotalLogicalProbability += P;
            }

            // Renormalize to fix leakage error if requested and valid
            if (renormalize && TotalLogicalProbability > 0.0)
            {
                const real_t InvTotal = 1.0 / TotalLogicalProbability;
                for (natural_t i = 0; i < LogicalDimSize; ++i)
                {
                    LogicalProbabilities[i] *= InvTotal;
                }
            }

            return LogicalProbabilities;
        }
    };
}