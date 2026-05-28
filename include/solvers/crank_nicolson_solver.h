#pragma once
#include "matrix_access.h"
#include "thomas_alogorithm.h"
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
    /// The Hamiltonian is provided explicitly as a matrix representation,
    /// typically originating from a spatial discretization of the kinetic
    /// and potential energy operators (e.g. finite-difference approximation
    /// of the Laplacian plus a potential term).
    ///
    /// Time integration is performed using the Crank–Nicolson scheme, which
    /// corresponds to an implicit midpoint discretization in time:
    ///
    ///   ( I + i·Δt/(2ℏ) · H ) · ψⁿ⁺¹ = ( I - i·Δt/(2ℏ) · H ) · ψⁿ
    ///
    /// This scheme is:
    ///  - second-order accurate in time,
    ///  - unconditionally stable,
    ///  - norm-preserving for Hermitian Hamiltonians,
    ///  - and therefore suitable for unitary quantum time evolution.
    ///
    /// In this implementation, the Hamiltonian is assumed to be time-independent
    /// (or piecewise constant in time). As a consequence, the Crank–Nicolson
    /// system matrices are constructed once and reused for each time step.
    ///
    /// If the Hamiltonian matrix is tridiagonal—as is the case for many
    /// one-dimensional finite-difference discretizations—the resulting linear
    /// system can be solved efficiently in O(N) time using the Thomas algorithm.

    /// @brief Callable object performing one Crank–Nicolson time step.
    ///
    /// @details
    /// The operator precomputes the Crank–Nicolson matrices A and B and
    /// applies them to advance a quantum state vector by a single time step.
    template<hilbert_space_t HilbertSpace,
        LinearSolverBackend Backend = LinearSolverBackend::ThomasTridiagonal>
    class CrankNicolsonSolver
    {
        static constexpr natural_t Dim = HilbertSpace::Dim;

        // Type alias for the matrix type used by the selected linear solver backend
        using solver_traits_t = LinearSolverTraits<Backend, Dim>;
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
        constexpr StateVector<HilbertSpace>
            operator()(const StateVector<HilbertSpace>& psi) const noexcept
        {
            // RHS = B · ψⁿ
            auto Rhs = multiply(m_B, psi);

            // Solve A · ψⁿ⁺¹ = RHS
            return LinearSolver<Backend, Dim>::
                template solve<HilbertSpace>(m_A, Rhs);
        }

        /// @brief  Construct/update the Crank–Nicolson system matrices A and B.
        /// @tparam Dim     Dimension of the Hilbert space.
        /// @param  hamiltonian  Hamiltonian operator of the system.
        /// @param  dt           Time step size.
        /// @param  A            Output matrix A = I + i·dt/(2ℏ)·H.
        /// @param  B            Output matrix B = I - i·dt/(2ℏ)·H.
        ///
        /// @details
        /// This function builds the two matrices required by the Crank–Nicolson
        /// time integration scheme. The matrices arise from the implicit midpoint
        /// discretization of the time-dependent Schrödinger equation.
        ///
        /// If the Hamiltonian matrix is tridiagonal, both A and B remain
        /// tridiagonal, enabling efficient O(N) time stepping.
        constexpr void updateMatrices(
            const matrix_type& H,
            real_t dt) noexcept
        {
            // i * dt / (2ℏ)
            // where we set ℏ = 1 in atomic units, so the factor simplifies to i * dt / 2
            const complex_t Factor(0.0, dt / 2.0);

            // Construct A = I + i·dt/(2ℏ)·H and B = I - i·dt/(2ℏ)·H
            for (natural_t i = 0; i < Dim; ++i)
            {
                for (natural_t j = 0; j < Dim; ++j)
                {
                    Matrix::set(m_A, i, j, 
                        Factor * Matrix::get(H, i, j));

                    Matrix::set(m_B, i, j, 
                        -Factor * Matrix::get(H, i, j));
                }

                // Add the identity matrix contribution to the main diagonals
                Matrix::set(m_A, i, i, 
                    complex_t::fromReal(1.0) + Matrix::get(m_A, i, i));
                Matrix::set(m_B, i, i, 
                    complex_t::fromReal(1.0) + Matrix::get(m_B, i, i));
            }
        }

    private:
        /// @brief  Helper function to compute the product of an instance of the helper tridiagonal matrix type
        ///         and a state vector.
        /// @tparam Dim     Dimension of the vector space.
        /// @param  M       Tridiagonal/square matrix.
        /// @param  psi    Input vector.
        /// @return         Resulting vector M · psi.
        constexpr StateVector<HilbertSpace> multiply(const matrix_type& M, const StateVector<HilbertSpace>& psi) const noexcept
        {
            StateVector<HilbertSpace> Result{ complex_t::zero() };

            if constexpr (Backend == LinearSolverBackend::ThomasTridiagonal)
            {
                Result[0] = M[MAINDIAGONAL][0] * psi[0] + M[SUPERDIAGONAL][0] * psi[1];

                for (natural_t i = 1; i < Dim - 1; ++i)
                {
                    Result[i] = M[SUBDIAGONAL][i] * psi[i - 1] +
                        M[MAINDIAGONAL][i] * psi[i] +
                        M[SUPERDIAGONAL][i] * psi[i + 1];
                }

                Result[Dim - 1] = M[SUBDIAGONAL][Dim - 1] * psi[Dim - 2] + M[MAINDIAGONAL][Dim - 1] * psi[Dim - 1];
            }
            else
            {
                for (natural_t i = 0; i < Dim; ++i) {
                    for (natural_t j = 0; j < Dim; ++j) {
                        Result[i] = Result[i] + M[i][j] * psi[j];
                    }
                }
            }
            return Result;
        }
    };
}