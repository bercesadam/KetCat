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
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
	struct EigenState
	{
		constexpr StateVector<HilbertSpace>
			operator()(unsigned int n) const noexcept
		{
			StateVector<HilbertSpace> Psi{};

			for (natural_t i = 0; i < HilbertSpace::Dim; ++i)
			{
				// Position (between Dirichlet boundaries)
				const real_t x = (i + 1) * HilbertSpace::dx;

				// Sin for the shape of the eigenstate
				const real_t Value = ConstexprMath::sin(
					n * ConstexprMath::Pi * x / HilbertSpace::Extent
				);

				Psi[i] = cplx_t(Value, 0.0);
			}

			Psi.normalize();
			return Psi;
		}
	};
}