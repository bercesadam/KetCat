#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"

namespace KetCat
{
	/// @brief Functor to generate the state vector corresponds with the Eigenstate
	///		   of the zero potential 1D box (with Dirichlet boundaries).
	/// @tparam Dim The dimension of the state vector to generate.
	/// @param n   Principal quantum number
	/// @param dx  Discretisation step
	/// @param L   Box length (w/ Dirichlet)
	template<natural_t Dim>
	struct EigenState
	{
		constexpr StateVector<InfiniteHilbertSpace1D<Dim>>
			operator()(unsigned int n, real_t dx, real_t L) const noexcept
		{
			StateVector<InfiniteHilbertSpace1D<Dim>> Psi{};

			for (natural_t i = 0; i < Dim; ++i)
			{
				// Position (between Dirichlet boundaries)
				const real_t x = (i + 1) * dx;

				// Sin for the shape of the eigenstate
				const real_t Value = ConstexprMath::sin(
					n * ConstexprMath::Pi * x / L
				);

				Psi[i] = cplx_t(Value, 0.0);
			}

			Psi.normalize(dx);
			return Psi;
		}
	};
}