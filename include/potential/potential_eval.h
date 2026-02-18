#pragma once
#include "core_types.h"

namespace KetCat
{
	/// @brief Concept to ensure that the PotentialFunctor is a callable
	///        object that takes a floating-point argument and returns a floating-point value.
	/// @tparam PotentialFunctor The type of the potential functor
	/// @tparam FloatType The floating-point type used for the input and output
	template <typename PotentialFunctor>
	concept potential_functor_t =
		std::is_class_v<PotentialFunctor> &&
		requires(PotentialFunctor obj, float_t val)
	{
		{ obj(val) } -> std::same_as<float_t>;
	};

	/// @brief Helper struct to evaluate multiple potential functors at a given position
	/// @details This struct allows for summing the contributions of multiple potential functors
	/// 	   at a specific position. It recursively evaluates each functor and accumulates their results.
	struct PotentialEvaluator
	{
		/// @brief Recursive operator() to evaluate and sum multiple potential functors at a given position
		/// @tparam First The first potential functor to evaluate
		/// @tparam Rest The remaining potential functors to evaluate
		template<potential_functor_t First, potential_functor_t... Rest>
		constexpr float_t operator()(float_t position, const First& first, const Rest&... rest) const
		{
			return first(position) + operator()(position, rest...);
		}

		/// @brief Base case for the recursive operator() when no more functors are left to evaluate
		constexpr float_t operator()(float_t position) const
		{
			return float_t{ 0 };
		}
	};
}
