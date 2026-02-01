#pragma once
#include "core_types.h"
#include "state_vector.h"

namespace KetCat
{
	/// @brief Functor to generate the state vector corresponds with the Eigenstate
	///		   of the zero potential 1D box (with Dirichlet boundaries).
	/// @tparam Dim The dimension of the state vector to generate.
	/// @param n   Principal quantum number
	/// @param dx  Discretisation step
	/// @param L   Box length (w/ Dirichlet)
	template<dimension_t Dim>
	struct EigenState
	{
		constexpr StateVector<Dim>
			operator()(unsigned int n, float_t dx, float_t L) const noexcept
		{
			StateVector<Dim> psi{};

			for (dimension_t i = 0; i < Dim; ++i)
			{
				// Position (between Dirichlet boundaries)
				const float_t x = (i + 1) * dx;

				// Sin for the shape of the eigenstate
				const float_t value = ConstexprMath::sin(
					n * ConstexprMath::Pi * x / L
				);

				psi[i] = cplx_t(value, 0.0);
			}

			psi.normalize_with_dx(dx);
			return psi;
		}
	};
}