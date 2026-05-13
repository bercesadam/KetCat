#pragma once
#include <functional>
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
        /// @brief Determine a basis state from a flat global state vector index, little-endian digit order.
        /// @details It's a mixed-radix (base-d) decoding to get the Qubit values at each position.
        ///          Example: For d=3, C=2, index 5 corresponds to basis state |2,1⟩ because 5 in base 3 is "21"
        /// @param globalIndex Global state vector index (0..d^C-1)
        /// @return Array of Qubit values (digits) corresponding to the basis state, ordered by position (little-endian)
        static constexpr qdit_list_t<QubitCount> decodeIndex(natural_t globalIndex) noexcept
        {
            qdit_list_t<QubitCount> BasisState{};
            for (natural_t i = 0; i < QubitCount; ++i)
            {
                BasisState[i] = globalIndex % LocalDim;
                globalIndex /= LocalDim;
            }
            return BasisState;
        }

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
        static constexpr bool isTargetsArrayValid(const qdit_list_t<K> targets) noexcept
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
        static constexpr bool isInTargets(natural_t position, const qdit_list_t<K>& targets) noexcept
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
        static constexpr natural_t nonTargetRank(natural_t position, const qdit_list_t<K>& targets) noexcept
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
                              const qdit_list_t<K>& targets) noexcept
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
        localTileIndex(const qdit_list_t<K>& targetDigits) noexcept
        {
            // little-endian
            natural_t idx = 0;
            natural_t mul = 1;
            for (natural_t i = 0; i < K; ++i)
            {
                idx += targetDigits[i] * mul;
                mul *= LocalDim;
            }
            return idx;
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
        static constexpr qdit_list_t<K> decodeLocalTileIndex(natural_t local) noexcept
        {
            std::array<natural_t, K> tdigits{};
            for (natural_t i = 0; i < K; ++i)
            {
                tdigits[i] = local % LocalDim;
                local /= LocalDim;
            }
            return tdigits;
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
        baseOffsetFromNonTargetIndex(natural_t nonTargetBasisIndex, const qdit_list_t<K>& targets) noexcept
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
        template <natural_t K>
        static constexpr void
        tileOperationsCore(const qdit_list_t<K> targetQdits,
                natural_t nonTargetBasisIndex,
                const std::function<void(natural_t, natural_t)>& operation) noexcept
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
            qdit_list_t<K> TileStrides{};
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
                const qdit_list_t<K>& targetQdits,
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
                    const qdit_list_t<K>& targetQdits,
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
        /// @brief Initialize all qdit registers based on a single-Qubit seed state |ψ⟩ ∈ ℂ^d,
        ///        creating the product state |Ψ⟩ = |ψ⟩^{⊗ C}.
        /// @param seed Single-Qubit state vector (size d)
        /// @return Product state vector for the entire register (size d^C)
        static constexpr StateVector<FullHilbertSpace>
            productStateFromSeed(const StateVector<OneQubitSpace>& seed) noexcept
        {
            StateVector<FullHilbertSpace> Result{};
            for (natural_t GlobalIndex = 0; GlobalIndex < FullDim; ++GlobalIndex)
            {
                auto Digits = decodeIndex(GlobalIndex);
                complex_t Amplitude = complex_t::fromReal(1.0);
                for (natural_t i = 0; i < QubitCount; ++i)
                {
                    Amplitude = Amplitude * seed[Digits[i]];
                }
                Result[GlobalIndex] = Amplitude;
            }
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
        template <natural_t K>
        static constexpr void applyHamiltonian(CrankNicolsonSolver<OperationSpace<K>>& solver,
                                            StateVector<FullHilbertSpace>& psi,
                                            qdit_list_t<K> targetQdits,
                                            const tridiagonal_matrix_t<OperationSpace<K>::Dim>& hamiltonian) noexcept
        {
            const natural_t BlockCount = blockCount<K>();

            StateVector<OperationSpace<K>> local{};
            StateVector<FullHilbertSpace> psiUpdated = psi;

            for (natural_t b = 0; b < BlockCount; ++b)
            {
                gatherTile<K>(psi, targetQdits, b, local);
                auto updatedLocal = solver(local);
                scatterTile<K>(psiUpdated, targetQdits, b, updatedLocal);
            }

            psi = psiUpdated;
        }

    private:
	    /// @brief  Compute the reduced density matrix of a single Qubit via partial trace.
        ///
        /// @param  psi         Global state vector (size d^C).
        /// @param  QubitIndex  Target Qubit whose local state is requested.
        ///
        /// @return             Reduced density matrix ρ_q = Tr_{¬q}(|ψ⟩⟨ψ|).
        ///
        /// @details
        /// The full pure-state density matrix is |ψ⟩⟨ψ|.  Tracing out all
        /// Qubits except Qubit q gives the d×d matrix:
        ///
        ///     ρ_q[a][b] = Σ_{env} ⟨a, env|ψ⟩⟨ψ|b, env⟩
        ///
        /// where the sum runs over all d^(C−1) environmental basis states.
        ///
        /// Equivalently, for each pair of global indices (i, j):
        /// - Decode both into Qubit digits.
        /// - Accept the pair only when all non-target digits agree
        ///   (partial-trace condition).
        /// - Accumulate:  ρ[dᵢ_q][dⱼ_q] += ψ[i] · ψ[j]*.
      /*  static constexpr DensityMatrix<LocalDim>
        reducedDensityMatrix(const StateVector<FullHilbertSpace>& psi,
                             natural_t                            QubitIndex) noexcept
        {
            DensityMatrix<LocalDim> Rho;
            Rho.setZero();

            for (natural_t i = 0; i < FullDim; ++i)
            {
                const auto Di = decodeIndex(i);

                for (natural_t j = 0; j < FullDim; ++j)
                {
                    const auto Dj = decodeIndex(j);

                    // Partial-trace condition: all non-target Qubits must match.
                    bool EnvMatch = true;
                    for (natural_t q = 0; q < QubitCount; ++q)
                    {
                        if (q == QubitIndex) continue;
                        if (Di[q] != Dj[q]) { EnvMatch = false; break; }
                    }
                    if (!EnvMatch) continue;

                    const natural_t A = Di[QubitIndex];
                    const natural_t B = Dj[QubitIndex];

                    Rho.m[A][B] += psi[i] * psi[j].conj();
                }
            }

            return Rho;
        }

        /// @brief  Extract a normalized ket from a rank-1 density matrix.
        ///
        /// @param  rho  Pure reduced density matrix (Tr(ρ²) ≈ 1).
        ///
        /// @return      Normalized state vector |ψ⟩ such that ρ ≈ |ψ⟩⟨ψ|.
        ///
        /// @details
        /// For a pure state ρ = |ψ⟩⟨ψ|, any non-zero column of ρ is
        /// proportional to |ψ⟩.  This function selects the column whose
        /// pivot row has the largest diagonal population (maximally stable
        /// numerical choice), then normalizes the result.
        ///
        /// Precondition: rho is positive-semidefinite with Tr(ρ²) ≈ 1.
        static constexpr StateVector<OneQubitSpace>
        pureStateVectorFromDensityMatrix(const DensityMatrix<LocalDim>& rho) noexcept
        {
            // Select pivot: row with largest diagonal element (highest population).
            natural_t Pivot   = 0;
            real_t    MaxPop  = rho.m[0][0].re;
            for (natural_t i = 1; i < LocalDim; ++i)
            {
                if (rho.m[i][i].re > MaxPop)
                {
                    MaxPop = rho.m[i][i].re;
                    Pivot  = i;
                }
            }

            // |ψ⟩ ∝ pivot-th column of ρ.
            StateVector<OneQubitSpace> Psi{};
            for (natural_t i = 0; i < LocalDim; ++i)
                Psi[i] = rho.m[i][Pivot];

            // Normalize.
            real_t Norm2 = real_t(0);
            for (natural_t i = 0; i < LocalDim; ++i)
                Norm2 += Psi[i].normSquared();

            if (Norm2 > real_t(0))
            {
                const real_t InvNorm = real_t(1) / ConstexprMath::sqrt(Norm2);
                for (natural_t i = 0; i < LocalDim; ++i)
                    Psi[i] = Psi[i] * InvNorm;
            }

            return Psi;
        }

    public:
        /// @brief  Extract the complete local-state description of a single Qubit.
        ///
        /// @param  psi         Global state vector representing the entire register.
        /// @param  QubitIndex  Index of the target Qubit (0 ≤ QubitIndex < C).
        ///
        /// @return             LocalStateInfo<d> containing:
        ///                     - The reduced density matrix ρ_q  (always valid).
        ///                     - The purity Tr(ρ²)               (always valid).
        ///                     - The entanglement classification  (Pure / Entangled).
        ///                     - A normalized ket |ψ⟩             (valid only when Pure).
        ///
        /// @details
        /// The function performs a partial trace over all Qubits except the
        /// target, producing the reduced density matrix ρ_q.
        ///
        /// The purity Tr(ρ²) distinguishes two physically distinct cases:
        ///
        ///  • Tr(ρ²) ≈ 1  →  Pure state.
        ///    The Qubit is unentangled from the rest of the register.
        ///    A unique (up to global phase) state vector |ψ⟩ exists and is
        ///    extracted via the dominant column of ρ.
        ///
        ///  • Tr(ρ²) < 1  →  Mixed / Entangled state.
        ///    The Qubit is entangled with the environment.
        ///    No single ket can describe its local state.
        ///    `pureStateVector` is left zeroed; use `rho` directly.
        ///
        /// The purity threshold is set to 1 − 1×10⁻¹⁰ to accommodate
        /// accumulated floating-point errors in constexpr arithmetic.
        ///
        /// Example usage:
        /// @code
        ///     auto Info = SubspaceHelper<3, 4>::extractLocalState(psi, 2);
        ///
        ///     if (Info.kind == LocalStateQualifier::Pure)
        ///     {
        ///         // Info.pureStateVector holds |ψ⟩ for Qubit 2.
        ///     }
        ///     else
        ///     {
        ///         // Qubit 2 is entangled — inspect Info.rho or Info.purityValue.
        ///     }
        /// @endcode
        static constexpr LocalStateInfo<LocalDim>
            extractLocalState(const StateVector<FullHilbertSpace>& psi,
                natural_t                            QubitIndex) noexcept
        {
            LocalStateInfo<LocalDim> Info{};

            if (QubitIndex >= QubitCount)
            {
                // Out-of-range: return zero-initialized descriptor.
                Info.rho.setZero();
                Info.purityValue = real_t(0);
                Info.kind = LocalStateQualifier::Entangled;
                return Info;
            }

            Info.rho = reducedDensityMatrix(psi, QubitIndex);
            Info.purityValue = purity(Info.rho);

            constexpr real_t PurityThreshold = real_t(1) - real_t(1e-10);

            if (Info.purityValue >= PurityThreshold)
            {
                Info.kind = LocalStateQualifier::Pure;
                Info.pureStateVector = pureStateVectorFromDensityMatrix(Info.rho);
            }
            else
            {
                Info.kind = LocalStateQualifier::Entangled;
                Info.pureStateVector = StateVector<OneQubitSpace>{};  // zeroed, intentionally invalid
            }

            return Info;
        }*/
    };

} 
