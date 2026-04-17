#pragma once
#include <array>
#include "constexprmath/constexpr_complex.h"
#include "constexprmath/constexpr_core_functions.h"
#include "constexprmath/constexpr_exp.h"


namespace KetCat
{
	/// @file
	/// @brief Core type aliases used across the project: sizes, number types, vectors and matrices.

	/// @brief Alias for unsigned integer types
	using natural_t = unsigned long;
	
	/// @brief Alias for floating-point type, enabling easy switching of precision if needed.
	using real_t = double;

	/// @brief Alias for complex numbers using the custom constexpr Complex type.
	using complex_t = ConstexprMath::Complex<real_t>;

	template<natural_t SpatialDimensions>
	using coordinate_t = std::array<natural_t, SpatialDimensions>;

	/// @brief State vector with compile-time fixed size (array of complex amplitudes).
	template<natural_t StateCount>
	using state_vector_t = std::array<complex_t, StateCount>;

	/// @brief Probability vector with compile-time fixed size (array of doubles).
	template<natural_t StateCount>
	using probability_vector_t = std::array<real_t, StateCount>;

	/// @brief Fixed-size list of qubit/qudit indices.
	/// @tparam QBitCountQDitCount  Number of qubits/qudits in the list.
	template<natural_t QBitCount>
	using qbit_list_t = std::array<natural_t, QBitCount>;
	template<natural_t QDitCount>
	using qdit_list_t = std::array<natural_t, QDitCount>;

	/// @brief Square matrix type 
	/// @tparam Rows  Number of rows and cols
	template<natural_t Dim>
	using matrix_t = std::array<std::array<complex_t, Dim>, Dim>;

	/// @brief Compact storage representation of a tridiagonal matrix for 1D Hamiltonians
	/// as we have useful information only in the three non-zero diagonals and this way we
	/// save memory by not storing Complex numbers of zero values
	/// The storage size of the upper and lower diagonial is the same as the main diagonal
	/// but the indicated values are never used.
	/// The major index 0 corresponds to the superdiagonal,
	///           index 1 corresponds to the main diagonal,
	///       and index 2 corresponds to the subdiagonal.
	template<natural_t Dim>
	using tridiagonal_matrix_t = std::array<std::array<complex_t, Dim>, 3U>;

	/// @brief Named constant indices for tridiagonal_matrix_t
	/// for convenience and intuitive usage.
	constexpr natural_t SUPERDIAGONAL = 0;
	constexpr natural_t MAINDIAGONAL = 1;
	constexpr natural_t SUBDIAGONAL = 2;
}

///@brief Tag struct to hold the spatial dimension's number
///       for the user defined literal.
struct DimensionTag
{
	KetCat::natural_t value;

	constexpr explicit DimensionTag(KetCat::natural_t v) : value(v) {}

	constexpr bool operator==(const DimensionTag& other) const noexcept
	{
		return value == other.value;
	}
};

///@brief User defined literal "_D" for elegant usage of spatial dimensions
constexpr DimensionTag operator ""_D(unsigned long long d)
{
	return DimensionTag{ static_cast<KetCat::natural_t>(d) };
}
