#pragma once
#include <array>
#include "constexprmath/constexpr_complex.h"
#include "constexprmath/constexpr_core_functions.h"

namespace KetCat
{
	/// @file
	/// @brief Core type aliases used across the project: sizes, complex number type, vectors and matrices.

	using dimension_t = std::size_t;
	using index_t = std::size_t;

	using float_t = double;

	using cplx_t = ConstexprMath::Complex<float_t>;

	/// @brief State vector with compile-time fixed size (array of complex amplitudes).
	template<dimension_t StateCount>
	using state_vector_t = std::array<cplx_t, StateCount>;

	/// @brief Probability vector with compile-time fixed size (array of doubles).
	template<dimension_t StateCount>
	using probability_vector_t = std::array<float_t, StateCount>;

	/// @brief Fixed-size list of qubit indices.
	/// @tparam QBitCount  Number of qubits in the list.
	template<index_t QBitCount>
	using qbit_list_t = std::array<dimension_t, QBitCount>;

	/// @brief Square matrix type 
	/// @tparam Rows  Number of rows and cols
	template<dimension_t Dim>
	using matrix_t = std::array<std::array<cplx_t, Dim>, Dim>;

	/// @brief Compact storage representation of a tridiagonal matrix for 1D Hamiltonians
	/// as we have useful information only in the three non-zero diagonals and this way we
	/// save memory by not storing Complex numbers of zero values
	/// The storage size of the upper and lower diagonial is the same as the main diagonal
	/// but the indicated values are never used.
	/// The major index 0 corresponds to the superdiagonal,
	///           index 1 corresponds to the main diagonal,
	///       and index 2 corresponds to the subdiagonal.
	template<dimension_t Dim>
	using tridiagonal_matrix_t = std::array<std::array<cplx_t, Dim>, 3U>;

	/// @brief Named constant indices for tridiagonal_matrix_t
	/// for convenience and intuitive usage.
	constexpr dimension_t SUPERDIAGONAL = 0;
	constexpr dimension_t MAINDIAGONAL = 1;
	constexpr dimension_t SUBDIAGONAL = 2;

}