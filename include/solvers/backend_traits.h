#pragma once
#include "core_types.h"


namespace KetCat
{
    /// @brief Enum for available linear solver backends
	/// to be used by the Crank–Nicolson solver
    enum class LinearSolverBackend
    {
        ThomasTridiagonal,
        GaussianElimination
    };


	/// @brief Primary template for linear solver traits, to be specialized by backend enumeration and dimension.
    template<LinearSolverBackend Backend, natural_t Dim>
    struct LinearSolverTraits;

    /// @brief Traits for Thomas tridiagonal solver.
    template<natural_t Dim>
    struct LinearSolverTraits<
        LinearSolverBackend::ThomasTridiagonal,
        Dim>
    {
        using matrix_type = tridiagonal_matrix_t<Dim>;
    };

    /// @brief Traits for dense Gaussian elimination solver.
    template<natural_t Dim>
    struct LinearSolverTraits<
        LinearSolverBackend::GaussianElimination,
        Dim>
    {
        using matrix_type = square_matrix_t<Dim>;
    };


	/// @brief Primary template for linear solvers, to be specialized by backend and dimension.
    template<LinearSolverBackend Backend, natural_t Dim>
    struct LinearSolver;
}
