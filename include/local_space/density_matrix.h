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
        static constexpr natural_t FullDim = ConstexprMath::pow(LevelCount, QubitCount);
        using FullHilbertSpace = FiniteHilbertSpace<FullDim>;

    public:
        /// @brief  Compute the reduced density matrix of a single Qubit via partial trace.
        ///
        /// @param  psi         Global state vector (size d^C).
        /// @param  qubitIndex  Target Qubit whose local state is requested.
        ///
        /// @return             Reduced density matrix ρ_q = Tr_{¬q}(|ψ⟩⟨ψ|).
        ///
        /// @details
        /// The full pure-state density matrix is |ψ⟩⟨ψ|. Tracing out all
        /// Qubits except Qubit q gives the d×d matrix:
        ///
        ///     ρ_q[a][b] = Σ_{env} ⟨a, env|ψ⟩⟨ψ|b, env⟩
        ///
        /// where the sum runs over all d^(C−1) environmental basis states.
        ///
        /// This implementation is optimized to O(d^(C+1)) complexity by directly 
        /// mapping the Little-Endian strides instead of performing a full O(d^2C) search.
        static constexpr Matrix<LevelCount>
            reducedDensityMatrix(const StateVector<FullHilbertSpace, QuantumPicture::Schrodinger>& psi,
                natural_t qubitIndex) noexcept
        {
            Matrix<LevelCount> Rho{};

            // Linear stride associated with the target qubit position in the global vector (Little-Endian)
            const natural_t QubitStride = ConstexprMath::pow(LevelCount, qubitIndex);

            // Total number of environmental configurations for the remaining (C - 1) qubits
            const natural_t EnvDim = ConstexprMath::pow(LevelCount, QubitCount - 1);

            // Iterate over all legetimate configurations of the environment
            for (natural_t env = 0; env < EnvDim; ++env)
            {
                // Embed the 'env' index into the global index, leaving the target qubit position as zero.
                // This follows the exact same mixed-radix logic as SubspaceHelper.
                natural_t BaseOffset = 0;
                natural_t Multiplier = 1;
                natural_t Rem = env;

                for (natural_t q = 0; q < QubitCount; ++q)
                {
                    if (q == qubitIndex)
                    {
                        Multiplier *= LevelCount; // Skip target qubit position
                    }
                    else
                    {
                        const natural_t Digit = Rem % LevelCount;
                        Rem /= LevelCount;
                        BaseOffset += Digit * Multiplier;
                        Multiplier *= LevelCount;
                    }
                }

                // With the environment fixed at BaseOffset, combine the local levels (A and B) of the target qubit
                for (natural_t A = 0; A < LevelCount; ++A)
                {
                    const natural_t GlobalI = BaseOffset + A * QubitStride;

                    for (natural_t B = 0; B < LevelCount; ++B)
                    {
                        const natural_t GlobalJ = BaseOffset + B * QubitStride;

                        // Accumulate: ρ[A][B] += ψ[GlobalI] · ψ[GlobalJ]*
                        Rho.at(A, B) = Rho.at(A, B) + psi[GlobalI] * psi[GlobalJ].conj();
                    }
                }
            }

            return Rho;
        }


        /// @brief  Calculates the purity of a density matrix: Tr(ρ²) ∈ (0, 1].
        ///
        /// @return      Tr(ρ²): equals 1 for a pure state, < 1 for a mixed (entangled) state.
        static constexpr real_t purity(Matrix<LevelCount> rho) noexcept
        {
            real_t P{};

            for (natural_t i = 0; i < LevelCount; ++i)
            {
                for (natural_t j = 0; j < LevelCount; ++j)
                {
                    // Tr(ρ²) = Σ_ij ρ_ij ρ_ji
                    P += (rho.at(i, j) * rho.at(j, i)).re;
                }
            }

            return P;
        }
    };
}
