#pragma once
#include "core_types.h"
#include "matrix_utils/matrix.h"


namespace KetCat
{
	/// @brief  Helper class to calculate the density matrix ρ of a given qubit from the global state vector
	///
	/// @details
	/// This class provides a static method to compute the reduced density matrix
	/// of a single Qubit by performing a partial trace over the rest of the system.
	/// which is used to extract local states and determine the entanglement structure of the system.
    /// 
	/// @tparam LevelCount  Local Hilbert-space dimension, corresponding
    ///         to the number of physical levels per qubit.
	/// @tparam QubitCount  Total number of qubits in the system
    template<natural_t LevelCount, natural_t QubitCount>
    class DensityMatrix
    {
        static constexpr natural_t FullDim = ConstexprMath::pow(LocalDim, QubitCount);
        using FullHilbertSpace = FiniteHilbertSpace<FullDim>;

    public:
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
        static constexpr Matrix<LevelCount>
            reducedDensityMatrix(const StateVector<FullHilbertSpace>& psi,
                natural_t qubitIndex) noexcept
        {
            Matrix<LevelCount> Rho{};

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
                        if (q == qubitIndex)
                        {
                            continue;
                        }

                        if (Di[q] != Dj[q])
                        {
                            EnvMatch = false;
                            break;
                        }
                    }

                    if (!EnvMatch)
                    {
                        continue;
                    }

                    const natural_t A = Di[QubitIndex];
                    const natural_t B = Dj[QubitIndex];

                    Rho.at(A, B) = Rho.at(A, B) + psi[i] * psi[j].conj();
                }
            }

            return Rho;
        }


    private:
        /// @brief Determine a basis state from a flat global state vector index, little-endian digit order.
        /// @details It's a mixed-radix (base-d) decoding to get the Qubit values at each position.
        ///          Example: For d=3, C=2, index 5 corresponds to basis state |2,1⟩ because 5 in base 3 is "21"
        /// @param globalIndex Global state vector index (0..d^C-1)
        /// @return Array of Qubit values (digits) corresponding to the basis state, ordered by position (little-endian)
        static constexpr qbit_list_t<QubitCount> decodeIndex(natural_t globalIndex) noexcept
        {
            qbit_list_t<QubitCount> BasisState{};
            for (natural_t i = 0; i < QubitCount; ++i)
            {
                BasisState[i] = globalIndex % LocalDim;
                globalIndex /= LocalDim;
            }
            return BasisState;
        }
    };
}
