#pragma once
#include "core_types.h"

namespace KetCat
{
    ///////// Linear Solver Backend

    /// @brief Enum for available linear solver backends
    /// to be used by the Crank–Nicolson solver.
    enum class LinearSolverBackend
    {
        ThomasTridiagonal,
        FiveBandGaussianElimination
    };


    ///////// Linear Solver Traits

    /// @brief Primary template for linear solver traits, to be specialized by backend enumeration and dimension.
    template<LinearSolverBackend Backend, natural_t LevelCount>
    struct LinearSolverTraits;

    /// @brief Traits for Thomas tridiagonal solver.
    template<natural_t LevelCount>
    struct LinearSolverTraits<LinearSolverBackend::ThomasTridiagonal, LevelCount>
    {
        using matrix_type = tridiagonal_matrix_t<LevelCount>;
    };

    /// @brief Traits for FiveBandGaussianElimination solver.
    template<natural_t LevelCount>
    struct LinearSolverTraits<LinearSolverBackend::FiveBandGaussianElimination, LevelCount>
    {
        using matrix_type = five_band_matrix_t<LevelCount>;
    };


    ///////// Linear Primary Template

    /// @brief Primary template for linear solvers, to be specialized by backend and dimension.
    template<LinearSolverBackend Backend, natural_t LevelCount>
    struct LinearSolver;


    ///////// System Size Traits

    /// @brief Primary template for system size traits, to be specialized by backend enumeration and dimension.
    template<natural_t LevelCount, LinearSolverBackend Backend>
    struct SystemSize;

    /// @brief Traits for Thomas Solver, where the dimension of the tridiagonal matrix
    /// (the length of the main diagonal) equals the number of levels in a basis set (one atom case)
    template<natural_t LevelCount>
    struct SystemSize<LevelCount, LinearSolverBackend::ThomasTridiagonal>
    {
        static constexpr natural_t value = LevelCount;
    };

    /// @brief Traits for the Five Band Matrix-based Gaussian Elimination vse, where the dimension of the matrix
    /// (the length of the main diagonal) equals to LevelCount², being the tensor product of two, one atom Hamiltonians
    template<natural_t LevelCount>
    struct SystemSize<LevelCount, LinearSolverBackend::FiveBandGaussianElimination>
    {
        static constexpr natural_t value = LevelCount * LevelCount;
    };
}
