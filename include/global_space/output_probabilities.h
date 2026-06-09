#pragma once
#include <bitset>
#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "local_space/neutral_atom_config.h"


namespace KetCat
{
    /// @brief  Extracts the measurement probabilities of the logical state vector.
    ///         Optionally renormalizes the results if probability leaked outside the computational subspace.
    /// @param  psi          The full global state vector (d^C dimensions).
    /// @param  logical0     Physical level index corresponding to logical |0>.
    /// @param  logical1     Physical level index corresponding to logical |1>.
    /// @param  renormalize  If true, rescales the probabilities to sum to 1.0 (conditional probability).
    /// @return              An array of size 2^QubitCount containing the probabilities in little-endian order.
    template <natural_t QubitCount, NeutralAtomTypeConfig Config>
    class LogicalProbabilityExtractor
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using GlobalStateManager = SubspaceHelper<ConfigType::LevelCount, QubitCount>;
        using FullStateVectorType = StateVector<typename GlobalStateManager::FullHilbertSpace>;
        
		// The size of the logical subspace is 2^QubitCount
        static constexpr natural_t LogicalDimSize = ConstexprMath::pow(natural_t(2), QubitCount);

    public:
		/// @brief Extracts the probabilities of the logical state vector from the full global state vector.
        static constexpr probability_vector_t<LogicalDimSize>
            extractLogicalProbabilities(const FullStateVectorType& psi, bool renormalize) noexcept
        {
            probability_vector_t<LogicalDimSize> LogicalProbabilities{};

            real_t TotalLogicalProbability = 0.0;

            // Step 1: Extract absolute probabilities and compute total subspace population
            for (natural_t logicalIdx = 0; logicalIdx < LogicalDimSize; ++logicalIdx)
            {
                std::bitset<QubitCount> bitstring(logicalIdx);

                natural_t GlobalIndex = 0;
                natural_t Multiplier = 1;

                for (natural_t i = 0; i < QubitCount; ++i)
                {
                    bool bit = bitstring[i];
                    natural_t PhysicalLevel = (bit ? ConfigType::Logical0Level : ConfigType::Logical1Level);

                    GlobalIndex += PhysicalLevel * Multiplier;
                    using ConfigType = std::remove_cvref_t<decltype(Config)>;
                    Multiplier *= ConfigType::LevelCount;
                }

                real_t P = psi[GlobalIndex].normSquared();
                LogicalProbabilities[logicalIdx] = P;
                TotalLogicalProbability += P;
            }

            // Step 2: Renormalize to fix leakage error if requested and valid
            if (renormalize && TotalLogicalProbability > 0.0)
            {
                real_t InvTotal = 1.0 / TotalLogicalProbability;
                for (natural_t i = 0; i < LogicalDimSize; ++i)
                {
                    LogicalProbabilities[i] *= InvTotal;
                }
            }

            return LogicalProbabilities;
        }
    };
}
