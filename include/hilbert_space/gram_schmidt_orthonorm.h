#pragma once
#include "basis_set.h"

namespace KetCat
{
    /// @brief Helper function, performs the good old BLAS-style axpy operation: Y = α·X + Y.
    /// @details Fundamental linear algebra building block used for both state vector 
    ///          manipulation and coefficient matrix updates. 
    /// @param alpha Scaling factor (complex).
    /// @param x     Input vector/array.
    /// @param y     Output vector/array to be modified in-place.
    /// @param size  Number of elements to process.
    template<natural_t Dim>
    constexpr void axpy(complex_t alpha, state_vector_t<Dim>& x, state_vector_t<Dim>& y) noexcept
    {
        for (natural_t k = 0; k < Dim; ++k)
        {
            y[k] = alpha * x[k] + y[k];
        }
    }

    /// @brief Gram-Schmidt Orthonormalization engine
    /// @tparam Dim Number of basis states to process.
    /// @tparam UseTwoPass If true, performs re-orthogonalization for better numerical stability.
    /// @details 
    /// This class calculates a transformation matrix L that maps a non-orthogonal 
    /// "raw" basis to an orthonormal set. It utilizes the Modified Gram-Schmidt (MGS) 
    /// algorithm, which is numerically superior to the Classical version.
	template<natural_t Dim, bool UseTwoPass = false>
    class Orthonormalizer
    {
        /// @brief Transformation recipe: orthonormal[i] = Σ L[i][j] * raw[j].
        matrix_t<Dim> m_Coefficients{};

        /// @brief Guard flag to ensure apply() is only called after training.
        bool m_isTrained = false;

    public:
        /// @brief Calculates and stores the orthonormalization recipe based on a sample basis.
        /// @tparam Space The underlying Hilbert space
        /// @param rawBasis The input basis set used to determine the transformation.
        /// @details
        /// Uses Modified Gram-Schmidt (MGS). By projecting the accumulating residual 
        /// (Current) instead of the raw input, we minimize the accumulation of 
        /// rounding errors.
        template<hilbert_space_t HilbertSpace>
        constexpr void learn(const BasisSet<HilbertSpace>& rawBasis)
        {
            auto StateVectors = rawBasis.getStateVectors();

            // Temporary buffer for orthonormalized vectors used during the j-loop projections.
            BasisSet<HilbertSpace> Temp;
            
            // Zero-initialize matrix and temporary basis
            m_Coefficients = {};
            Temp = {};

            // Outer loop thru all basis states
            for (natural_t i = 0; i < Dim; ++i)
            {
                auto Current = StateVectors[i]; 
                m_Coefficients[i][i] = 1.0; 

                // --- Modified Gram-Schmidt Step --- 
                for (natural_t j = 0; j < i; ++j)
                {
                    /// Pass 1: Project the residual 'Current' onto the already normalized Temp[j].
                    // MGS Improvement: Using 'Current' instead of 'rawBasis[i].m_Psi' 
                    // significantly improves numerical stability.
                    complex_t innerProduct1 = Temp[j].innerProduct(Current);

                    axpy(-innerProduct1, m_Coefficients[j], m_Coefficients[i]);
                    axpy(-innerProduct1, Temp[j], Current);

                    // Optional Pass 2: Re-orthogonalization ("Twice is Enough").
                    // Effectively cleans up remaining machine-precision leakage.
                    if constexpr (UseTwoPass)
                    {
                        complex_t innerProduct2 = Temp[j].innerProduct(Current);
                        axpy(-innerProduct2, m_Coefficients[j], m_Coefficients[i]);
                        axpy(-innerProduct2, Temp[j], Current);
                    }
                }

                // --- Normalization Step ---
                real_t Norm = ConstexprMath::sqrt(Current.innerProduct(Current).re);

                // Safety measure, handle linear dependency
                if (Norm < 1e-12)
                {
                    [] { throw "Linear dependence in basis set!"; }();
                }
                
                real_t InvNorm = 1.0 / Norm;
                
                // Scale the coefficients for the final recipe.
                for (natural_t k = 0; k <= i; ++k)
                {
                    m_Coefficients[i][k] =  m_Coefficients[i][k] * InvNorm;
                }
                
                // Scale the wavefunction for future j-loop projections.
                for (natural_t k = 0; k < Dim; ++k)
                {
                    Current[k] = Current[k] * InvNorm;

                }

                Temp[i] = Current;
            }

            m_isTrained = true;
        }

        /// @brief Applies the learned transformation recipe to a new basis set.
        /// @tparam Space Target Hilbert space.
        /// @param rawBasis A new basis set (must follow the same ordering as learn()).
        /// @return A transformed, orthonormalized BasisSet.
        template<hilbert_space_t HilbertSpace>
        constexpr BasisSet<HilbertSpace> apply(const const BasisSet<HilbertSpace>& rawBasis) const
        {
            if (!m_isTrained) 
            {
                [] { throw "Orthonormalizer: learn() must be called before apply()!"; }();
            }

            BasisSet<HilbertSpace> Transformed;

            for (size_t i = 0; i < N; ++i)
            {
                Transformed[i].m_Psi = {};
                Transformed[i].m_Energy = rawBasis[i].m_Energy;

                for (size_t j = 0; j <= i; ++j)
                {
                    axpy(m_Coefficients[i][j], rawBasis[j], Transformed[i].m_Psi)
                }
            }

            return Transformed;
        }
    };
}


