#pragma once
#include <bitset>

#include "hilbert_space/hilbert.h"        
#include "hilbert_space/state_vector.h" 
#include "solvers/crank_nicolson_solver.h"


namespace KetCat
{
    template<natural_t _LocalQditDim, natural_t _QubitCount>
    class SubspaceHelper
    {
        static_assert(_LocalQditDim >= 2, "Local dimension must be >= 2");
        static_assert(_QubitCount >= 1, "Qubit count must be >= 1");

    public:
        static constexpr natural_t LocalDim = _LocalQditDim;
        static constexpr natural_t QubitCount = _QubitCount;
        static constexpr natural_t FullDim = ConstexprMath::pow(LocalDim, QubitCount);

        using FullHilbertSpace = FiniteHilbertSpace<FullDim>;
        using OneQubitSpace   = FiniteHilbertSpace<LocalDim>;

        template<natural_t TargetQdits>
        using OperationSpace = FiniteHilbertSpace<ConstexprMath::pow(LocalDim, TargetQdits)>;

    private:
        /// @brief  Number of tiles produced when selecting K target Qubits.
        ///
        /// @tparam K  Number of Qubits that form the local tile.
        ///
        /// @return Number of distinct tiles in the global state vector.
        ///
        /// @details
        /// When K Qubits are selected as targets, the remaining (C − K) Qubits
        /// determine which slice of the global state vector is accessed.
        /// Each configuration of the non-target Qubits corresponds to one tile.
        ///
        /// The number of such configurations is:
        ///
        ///     d^(C − K)
        ///
        /// where:
        /// - C is the total number of Qubits
        /// - d is the local dimension of each Qubit.
        ///
        /// Each tile contains exactly d^K amplitudes.
        template <natural_t K>
        static constexpr natural_t blockCount() noexcept
        {
            static_assert(K <= QubitCount, "K cannot exceed the number of Qubits");
            return ConstexprMath::pow(LocalDim, QubitCount - K);
        }

        /// @brief  Validate the list of target Qubit indices.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  targets  Array containing the indices of the selected Qubits.
        ///
        /// @return true if all indices are valid and unique.
        ///
        /// @details
        /// This function verifies two conditions:
        ///
        /// 1. All indices lie within the valid Qubit range [0, C).
        /// 2. No index appears more than once.
        ///
        /// The function performs a simple O(K²) uniqueness check, which is
        /// acceptable because K is typically small (e.g. 1–3 in most
        /// quantum gate operations).
        ///
        /// This validation is primarily intended for compile-time or
        /// initialization checks when defining a subsystem.
        template <natural_t K>
        static constexpr bool isTargetsArrayValid(const qbit_list_t<K> targets) noexcept
        {
            // Check that all indices are within the valid range [0, C)
            for (natural_t i = 0; i < K; ++i)
            {
                if (targets[i] >= QubitCount)
                {
                    return false;
                }
            }

            // Check for uniqueness of indices
            for (natural_t i = 0; i < K; ++i)
            {
                for (natural_t j = i + 1; j < K; ++j)
                {
                    if (targets[i] == targets[j])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        // @brief  Check whether a global Qubit position belongs to the target set.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  pos      Global Qubit position in the range [0, C).
        /// @param  targets  Array of target Qubit indices.
        ///
        /// @return true if pos is one of the target Qubits.
        ///
        /// @details
        /// This helper function performs a linear membership test over the
        /// target index array. It is used by several algorithms that need
        /// to distinguish between:
        ///
        /// - target Qubits (which vary inside the tile)
        /// - non-target Qubits (which select the tile).
        template <natural_t K>
        static constexpr bool isInTargets(natural_t position, const qbit_list_t<K>& targets) noexcept
        {
            for (natural_t i = 0; i < K; ++i)
            {
                if (targets[i] == position)
                {
                    return true;
                }
            }
            return false;
        }

        /// @brief  Rank of a non-target Qubit position among all non-target Qubits.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  pos      Global Qubit position whose rank is requested.
        /// @param  targets  Array of target Qubit indices.
        ///
        /// @return Zero-based rank of the non-target position.
        ///
        /// @details
        /// When scanning the global Qubit positions in ascending order
        /// (0 … C−1), this function counts how many *non-target* Qubits
        /// appear before the given position.
        ///
        /// The resulting rank determines which digit of the
        /// `nonTargetBasisIndex` corresponds to this Qubit.
        ///
        /// Example:
        ///
        ///     positions:   0 1 2 3 4
        ///     targets:       ^   ^
        ///
        /// non-target positions: {0,2,4}
        ///
        /// ranks:
        ///
        ///     pos 0 → rank 0
        ///     pos 2 → rank 1
        ///     pos 4 → rank 2
        template <natural_t K>
        static constexpr natural_t nonTargetRank(natural_t position, const qbit_list_t<K>& targets) noexcept
        {
            natural_t Rank = 0;
            for (natural_t i = 0; i < position; ++i)
            {
                if (!isInTargets(i, targets))
                {
                    ++Rank;
                }
            }
            return Rank;
        }

        /// @brief  Extract the base-d digit for a specific non-target Qubit.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  nonTargetBasisIndex Index describing the configuration of non-target Qubits.
        /// @param  pos      Global Qubit position whose digit is requested.
        /// @param  targets  Array of target Qubit indices.
        ///
        /// @return Base-d digit corresponding to Qubit position `pos`.
        ///
        /// @details
        /// The nonTargetBasisIndex encodes the basis state of
        /// the non-target Qubits using base-d digits.
        ///
        /// This function determines which digit of that number corresponds
        /// to the given global Qubit position.
        ///
        /// Precondition:
        ///
        ///     pos ∉ targets
        ///
        /// The digit is obtained by:
        ///
        /// 1. Determining the rank of the non-target Qubit.
        /// 2. Extracting the corresponding base-d digit from the index.
        template <natural_t K>
        static constexpr natural_t
        digitAtPosFromBlockId(natural_t nonTargetBasisIndex, natural_t position,
                              const qbit_list_t<K>& targets) noexcept
        {
            const natural_t Rank = nonTargetRank(position, targets);
            return (nonTargetBasisIndex / ConstexprMath::pow(LocalDim, Rank)) % LocalDim;
        }

        /// @brief  Encode target-Qubit digits into a local tile index.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  tdigits  Base-d digits representing the states of the
        ///                  target Qubits in the order specified by `targets`.
        ///
        /// @return Linear index inside the tile.
        ///
        /// @details
        /// The digits are interpreted in little-endian order with respect
        /// to the provided target ordering:
        ///
        ///     index = t0 + t1·d + t2·d² + ...
        ///
        /// where:
        ///
        /// - ti is the state of the i-th target Qubit
        /// - d is the local Qubit dimension.
        ///
        /// This mapping converts the K-dimensional target subsystem
        /// into a contiguous one-dimensional tile index.
        template <natural_t K>
        static constexpr natural_t
        localTileIndex(const qbit_list_t<K>& targetDigits) noexcept
        {
            // little-endian
            natural_t Index = 0;
            natural_t Multiplier = 1;
            for (natural_t i = 0; i < K; ++i)
            {
                Index += targetDigits[i] * Multiplier;
                Multiplier *= LocalDim;
            }
            return Index;
        }

        /// @brief  Decode a local tile index into base-d digits.
        ///
        /// @tparam K     Number of target Qubits.
        /// @param  local Linear index within the tile.
        ///
        /// @return Array containing the base-d digits for each target Qubit.
        ///
        /// @details
        /// This is the inverse operation of `localTileIndex`.
        ///
        /// The function decomposes the linear tile index into K base-d digits
        /// representing the states of the target Qubits:
        ///
        ///     local = t0 + t1·d + t2·d² + ...
        ///
        /// The digits are returned in little-endian order matching the
        /// ordering of the `targets` array.
        template <natural_t K>
        static constexpr qbit_list_t<K> decodeLocalTileIndex(natural_t local) noexcept
        {
            std::array<natural_t, K> Digits{};
            for (natural_t i = 0; i < K; ++i)
            {
                Digits[i] = local % LocalDim;
                local /= LocalDim;
            }
            return Digits;
        }

        /// @brief  Compute the base global index of a tile.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  nonTargetBasisIndex  Basis index describing the configuration
        ///                              of the non-target Qubits.
        /// @param  targets  Array of target Qubit indices.
        ///
        /// @return Base offset of the tile within the global state vector.
        ///
        /// @details
        /// The global basis index is constructed by embedding the base-d
        /// digits of `nonTargetBasisIndex` into the non-target Qubit positions
        /// of the full system.
        ///
        /// Target Qubit positions are left as zeros because their values
        /// will be added later when iterating inside the tile.
        ///
        /// Conceptually this function constructs:
        ///
        ///     | q(C−1) ... q0 >
        ///
        /// where:
        ///
        /// - target Qubits are set to 0
        /// - non-target Qubits are taken from `nonTargetBasisIndex`.
        ///
        /// The result is the starting index of the tile in the global state vector.
        template <natural_t K>
        static constexpr natural_t
        baseOffsetFromNonTargetIndex(natural_t nonTargetBasisIndex, const qbit_list_t<K>& targets) noexcept
        {
            natural_t Offset = 0;
            natural_t Multiplier = 1;
            natural_t rem = nonTargetBasisIndex;
            for (natural_t pos = 0; pos < QubitCount; ++pos)
            {
                if (isInTargets(pos, targets))
                {
                    Multiplier *= LocalDim;
                }
                else
                {
                    const natural_t Digit = rem % LocalDim;
                    rem /= LocalDim;
                    Offset += Digit * Multiplier;
                    Multiplier *= LocalDim;
                }
            }
            return Offset;
        }

        /// @brief  Core traversal routine for a target-Qubit tile.
        ///
        /// @tparam K  Number of target Qubits.
        /// @param  targetQdits           Indices of the target Qubits.
        /// @param  nonTargetBasisIndex   Index selecting the tile.
        /// @param  operation             Callback invoked for each tile element.
        ///
        /// @details
        /// This function performs the common traversal used by both
        /// tile gathering and scattering operations.
        ///
        /// The algorithm:
        ///
        /// 1. Computes the base global offset of the tile using
        ///    `nonTargetBasisIndex`.
        /// 2. Iterates over all d^K basis states of the target Qubits.
        /// 3. Constructs the corresponding global index.
        /// 4. Calls the provided operation with:
        ///
        ///        (localTileIndex, globalIndex)
        ///
        /// This abstraction allows different operations (copy, transform,
        /// accumulate, etc.) to reuse the same indexing logic.
        template <natural_t K, typename Op>
        static constexpr void
            tileOperationsCore(
                const qbit_list_t<K>& targetQdits,
                natural_t nonTargetBasisIndex,
                Op&& operation) noexcept
        {
            // `targetQdits` is a function parameter and may not be a compile-time
            // constant. Use a runtime check instead of static_assert which
            // requires a constant expression.
            if (!isTargetsArrayValid(targetQdits))
            {
                return; // invalid targets: nothing to do
            }

            // Compute the linear strides associated with a Qubit position.
            // The global basis index of a multi-Qubit state is represented as a base-d number:
            // index = q0 + q1·d + q2·d² + ... + q(C−1)·d^(C−1)
            // This value represents how much the global linear index changes when the state of Qubit `i` is increased by one.
            qbit_list_t<K> TileStrides{};
            for (natural_t i = 0; i < K; ++i)
            {

                TileStrides[i] = ConstexprMath::pow(LocalDim, targetQdits[i]);
            }

            //  Size of a local tile corresponding to K target Qubits.
            const natural_t TileSize = ConstexprMath::pow(LocalDim, K);
            // Base global index of the tile, determined by the configuration of non-target Qubits.
            const natural_t BaseOffset = baseOffsetFromNonTargetIndex(nonTargetBasisIndex, targetQdits);

            // Iterate over all local indices of the tile and compute the corresponding global index.
            // Then apply the provided operation for each pair of local and global indices.
            for (natural_t LocalIndex = 0; LocalIndex < TileSize; ++LocalIndex)
            {
                const auto TileDigits = decodeLocalTileIndex<K>(LocalIndex);
                natural_t GlobalIndex = BaseOffset;
                for (natural_t i = 0; i < K; ++i)
                {
                    GlobalIndex += TileDigits[i] * TileStrides[i];
                }

                operation(LocalIndex, GlobalIndex);
            }
        }

        /// @brief  Gather a tile of amplitudes from the global state vector.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  fullSpace     Global state vector.
        /// @param  targetQdits  Target Qubit indices.
        /// @param  nonTargetBasisIndex Index selecting the tile.
        /// @param  out      Output tile of size d^K.
        ///
        /// @details
        /// This function extracts the amplitudes corresponding to the
        /// selected target Qubits while the remaining Qubits are fixed
        /// according to `block_id`.
        ///
        /// The result is a contiguous tile containing all d^K amplitudes
        /// of the target subsystem.
        template <natural_t K>
        static constexpr void
        gatherTile(const StateVector<FullHilbertSpace>& fullSpace,
                const qbit_list_t<K>& targetQdits,
                natural_t nonTargetBasisIndex,
                StateVector<FiniteHilbertSpace<ConstexprMath::pow(LocalDim, K)>>& out) noexcept
        {
            tileOperationsCore<K>(targetQdits, nonTargetBasisIndex,
                [&](natural_t LocalIndex, natural_t GlobalIndex)
                {
                    out[LocalIndex] = fullSpace[GlobalIndex];
                });
        }

        /// @brief  Scatter a tile of amplitudes back into the global state vector.
        ///
        /// @tparam K        Number of target Qubits.
        /// @param  fullSpace     Global state vector.
        /// @param  targetQdits  Target Qubit indices.
        /// @param  nonTargetBasisIndex Index selecting the tile.
        /// @param  in       Tile containing d^K amplitudes.
        ///
        /// @details
        /// This function performs the inverse operation of `gatherTile`.
        /// The provided tile is written back into the global state vector
        /// at the positions corresponding to the specified tile.
        ///
        /// Each amplitude of the tile is mapped to its corresponding
        /// global basis index.
        template <natural_t K>
        static constexpr void
        scatterTile(StateVector<FullHilbertSpace>& fullSpace,
                    const qbit_list_t<K>& targetQdits,
                    natural_t nonTargetBasisIndex,
                    const StateVector<FiniteHilbertSpace<ConstexprMath::pow(LocalDim, K)>>& in) noexcept
        {
            tileOperationsCore<K>(targetQdits, nonTargetBasisIndex,
                [&](natural_t LocalIndex, natural_t GlobalIndex)
                {
                    fullSpace[GlobalIndex] = in[LocalIndex];
                });
        }

    public:
		/// @brief Initialize the global state vector to a specific computational basis state defined by a bitstring.
		/// @param bitstring   The bitstring representing the desired basis state, where each bit corresponds to a Qubit (0 or 1).
		/// @param logical0    Physical level corresponding to the logical '0' state
        /// @param logical1    Physical level corresponding to the logical '1' state
		/// @return            A state vector of size d^C with amplitude 1 at the index corresponding
        ///                    to the specified bitstring and 0 elsewhere.
        static constexpr StateVector<FullHilbertSpace>
            basisStateFromBitstring(std::bitset<QubitCount> bitstring,
                natural_t logical0, natural_t logical1) noexcept
        {
            StateVector<FullHilbertSpace> Result{};

            natural_t GlobalIndex = 0;
            natural_t Multiplier = 1;

            for (natural_t i = 0; i < QubitCount; ++i)
            {
                bool bit = bitstring[i];
                natural_t PhysicalLevel = (bit ? logical1 : logical0);

                GlobalIndex += PhysicalLevel * Multiplier;
                Multiplier *= LocalDim;
            }   

			// Fill the global state vector with zeros except for the specified basis state which is set to amplitude 1.
            Result[GlobalIndex] = complex_t::fromReal(1.0);

            return Result;
        }

        /// @brief  Apply a K-Qubit operation defined by a tridiagonal Hamiltonian
        ///         to the specified target Qubits in the global state vector.
        /// @tparam K        Number of target Qubits.
        /// @param  psi             Global state vector representing the entire register.
        /// @param  targetQdits     Indices of the target Qubits to which the operation is applied.
        /// @param  hamiltonian     Tridiagonal matrix defining the K-Qubit operation in the local tile space.
        /// @param  dt              Time step size for the Crank–Nicolson evolution.
        /// @return                 Updated global state vector after applying the operation.
        /// @details
        /// This function performs the following steps:
        /// 1. Iterates over all tiles corresponding to the selected target qubits.
        /// 2. For each tile:
        ///    a. Gathers the relevant amplitudes from the global state vector into a local tile vector.
        ///    b. Applies the Crank–Nicolson time evolution using the provided Hamiltonian.
        ///    c. Scatters the updated tile amplitudes back into the global state vector.
        template <natural_t K, LinearSolverBackend L>
        static constexpr void performTimeEvolution(CrankNicolsonSolver<LocalDim, L>& solver,
                                            StateVector<FullHilbertSpace>& psi,
                                            qbit_list_t<K> targetQdits) noexcept
        {
            const natural_t BlockCount = blockCount<K>();

            StateVector<FullHilbertSpace> psiUpdated = psi;

            for (int b = 0; b < BlockCount; ++b)
            {
                StateVector<OperationSpace<K>> local{};
                gatherTile<K>(psi, targetQdits, b, local);
                auto updatedLocal = solver(local);
                scatterTile<K>(psiUpdated, targetQdits, b, updatedLocal);
            }

            psi = psiUpdated;
        }
    };
} 
