#pragma once
#include "matrix_access.h"
#include "thomas_algorithm.h"
#include "gauss_elimination.h"


namespace KetCat
{
    /// @brief Crank–Nicolson solver for the time-dependent Schrödinger equation.
    ///
    /// @details
    /// This file implements a numerical time evolution scheme for quantum
    /// state vectors based on the Crank–Nicolson method. The solver advances
    /// a discretized wavefunction by solving the time-dependent Schrödinger equation
    ///
    ///   iℏ ∂ψ(t)/∂t = H ψ(t)
    ///
    /// where ψ(t) is a state vector in a finite-dimensional Hilbert space and
    /// H is the Hamiltonian operator of the system.
    ///
    /// The Hamiltonian is provided explicitly as a matrix representation. Depending
    /// on the chosen backend, it targets specific grid-based structures:
    ///  - ThomasTridiagonal: Designed for 1-dimensional finite-difference stencils.
    ///  - FiveBandElimination: Designed for 2-dimensional Hamiltonians originating from the 
    ///    2D Laplace operator, which can be expressed as the tensor product of two 
    ///    independent 1D tridiagonal sub-operators.
    ///
    /// Time integration is performed using the Crank–Nicolson scheme, which
    /// corresponds to an implicit midpoint discretization in time:
    ///
    ///   ( I + i·Δt/(2ℏ) · H ) · ψⁿ⁺¹ = ( I - i·Δt/(2ℏ) · H ) · ψⁿ
    ///
    /// This scheme is second-order accurate in time, unconditionally stable,
    /// and norm-preserving for Hermitian Hamiltonians.
    ///
    /// @tparam HilbertSpace The underlying quantum Hilbert space structure.
    /// @tparam Backend      The linear solver backend to utilize (Tridiagonal or FiveBandElimination).
    template<natural_t LevelCount, LinearSolverBackend Backend>
    class CrankNicolsonSolver
    {
        static constexpr natural_t SystemLevelCount = SystemSize<LevelCount, Backend>::value;

        // Type alias for the matrix type used by the selected linear solver backend
        using solver_traits_t = LinearSolverTraits<Backend, LevelCount>;
        using matrix_type = typename solver_traits_t::matrix_type;

        // Alias for matrix access traits
        using Matrix = MatrixAccess<matrix_type>;

        // Precomputed matrices
        matrix_type m_A;
        matrix_type m_B;

    public:
        /// @brief  Advances the state vector by one time step.
        /// @param  psi     State vector at time step n.
        /// @return         State vector at time step n+1.
        ///
        /// @details
        /// The function computes the right-hand side B · ψⁿ and then solves
        /// the linear system A · ψⁿ⁺¹ = RHS, resulting in unitary time evolution
        /// using the selected linear solver backend.
        template<hilbert_space_t HilbertSpace>
        constexpr StateVector<HilbertSpace>
            operator()(const StateVector<HilbertSpace>& psi) const noexcept
        {
            // RHS = B · ψⁿ
            auto Rhs = multiply(m_B, psi);

            // Solve A · ψⁿ⁺¹ = RHS
            return LinearSolver<Backend, LevelCount>::
                template solve<HilbertSpace>(m_A, Rhs);
        }

        /// @brief  Construct/update the Crank–Nicolson system matrices A and B.
        /// @param  H            Hamiltonian operator of the system.
        /// @param  dt           Time step size.
        ///
        /// @details
        /// This function builds the two matrices required by the Crank–Nicolson
        /// time integration scheme. The loops are highly optimized via constexpr branch
        /// selection to ensure only the occupied bands are mapped from the source Hamiltonian.
        constexpr void updateMatrices(const matrix_type& H, real_t dt) noexcept
        {
            const complex_t Factor(0.0, dt / 2.0);

            // Construct A = I + i·dt/(2ℏ)·H and B = I - i·dt/(2ℏ)·H with
            // optimized O(N) five band matrix assembly loop targeting only populated diagonals

            if constexpr (Backend == LinearSolverBackend::ThomasTridiagonal)
            {
                for (natural_t i = 0; i < SystemLevelCount; ++i)
                {
                    // Main diagonal processing including identity addition
                    Matrix::set(m_A, i, i, complex_t::fromReal(1.0) + Factor * Matrix::get(H, i, i));
                    Matrix::set(m_B, i, i, complex_t::fromReal(1.0) - Factor * Matrix::get(H, i, i));

                    // Immediate sub and super diagonals
                    if (i > 0)
                    {
                        Matrix::set(m_A, i, i - 1, Factor * Matrix::get(H, i, i - 1));
                        Matrix::set(m_B, i, i - 1, -Factor * Matrix::get(H, i, i - 1));
                    }
                    if (i < SystemLevelCount - 1)
                    {
                        Matrix::set(m_A, i, i + 1, Factor * Matrix::get(H, i, i + 1));
                        Matrix::set(m_B, i, i + 1, -Factor * Matrix::get(H, i, i + 1));
                    }
                }
            }
            else if constexpr (Backend == LinearSolverBackend::FiveBandGaussianElimination)
            {
                
                for (natural_t i = 0; i < SystemLevelCount; ++i)
                {
                    // Main diagonal processing including identity addition
                    Matrix::set(m_A, i, i,
                        complex_t::fromReal(1.0) + Factor * Matrix::get(H, i, i));
                    Matrix::set(m_B, i, i,
                        complex_t::fromReal(1.0) - Factor * Matrix::get(H, i, i));

                    // Immediate sub and super diagonals
                    if (i > 0)
                    {
                        Matrix::set(m_A, i, i - 1,
                            Factor * Matrix::get(H, i, i - 1));
                        Matrix::set(m_B, i, i - 1,
                           -Factor * Matrix::get(H, i, i - 1));
                    }
                    if (i < SystemLevelCount - 1)
                    {
                        Matrix::set(m_A, i, i + 1,
                            Factor * Matrix::get(H, i, i + 1));
                        Matrix::set(m_B, i, i + 1,
                           -Factor * Matrix::get(H, i, i + 1));
                    }

                    // Far upper and lower diagonals
                    if (i >= LevelCount) 
                    {
                        Matrix::set(m_A, i, i - LevelCount,
                            Factor * Matrix::get(H, i, i - LevelCount)); // i - 2 helyett
                        Matrix::set(m_B, i, i - LevelCount,
                           -Factor * Matrix::get(H, i, i - LevelCount));
                    }
                    if (i < SystemLevelCount - LevelCount)
                    {
                        Matrix::set(m_A, i, i + LevelCount,
                            Factor * Matrix::get(H, i, i + LevelCount)); // i + 2 helyett
                        Matrix::set(m_B, i, i + LevelCount,
                           -Factor * Matrix::get(H, i, i + LevelCount));
                    }
                }
            }
        }


    private:
        /// @brief  Computes the matrix-vector product of a banded matrix and a state vector.
        /// @param  M       Banded matrix storage framework instance.
        /// @param  psi     Input quantum state vector.
        /// @return         Resulting state vector M · psi.
        ///
        /// @details
        /// Dispatches at compile time to the most efficient multiplication algorithm.
        /// Eradicates dense O(N^2) loops in favor of strict O(N) banded traversals.
        template<hilbert_space_t HilbertSpace>
        constexpr StateVector<HilbertSpace>
            multiply(const matrix_type& M,
                const StateVector<HilbertSpace>& psi) const noexcept
        {
            StateVector<HilbertSpace> Result{ complex_t::zero() };

            if constexpr (Backend == LinearSolverBackend::ThomasTridiagonal)
            {
                // Unrolled O(N) tridiagonal matrix-vector multiplication skipping outer zero fields
                for (natural_t i = 0; i < LevelCount; ++i)
                {
                    Result[i] = M[MAINDIAGONAL][i] * psi[i];

                    if (i > 0)
                    {
                        Result[i] = Result[i] + M[SUBDIAGONAL][i] * psi[i - 1];
                    }

                    if (i < LevelCount - 1)
                    {
                        Result[i] = Result[i] + M[SUPERDIAGONAL][i] * psi[i + 1];
                    }
                }
            }
            else if constexpr (Backend == LinearSolverBackend::FiveBandGaussianElimination)
            {
                // Unrolled O(N) pentadiagonal matrix-vector multiplication skipping outer zero fields
                for (natural_t i = 0; i < SystemLevelCount; ++i) {
                    Result[i] = M[MAINDIAGONAL][i] * psi[i];

                    // Közvetlen alsó szomszéd (y-irányú átmenetek)
                    if (i > 0) {
                        Result[i] = Result[i] + M[SUBDIAGONAL][i] * psi[i - 1];
                    }
                    // Távoli alsó sáv (x-irányú átmenetek, távolság: LevelCount)
                    if (i >= LevelCount) {
                        Result[i] = Result[i] + M[LOWER_FAR][i] * psi[i - LevelCount];
                    }
                    // Közvetlen felső szomszéd (y-irányú átmenetek)
                    if (i < SystemLevelCount - 1) {
                        Result[i] = Result[i] + M[SUPERDIAGONAL][i] * psi[i + 1];
                    }
                    // Távoli felső sáv (x-irányú átmenetek, távolság: LevelCount)
                    if (i < SystemLevelCount - LevelCount) {
                        Result[i] = Result[i] + M[UPPER_FAR][i] * psi[i + LevelCount];
                    }
                }
            }

            return Result;
        }
    };
}
