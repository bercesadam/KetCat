#pragma once
#include "core_types.h"
#include "state_vector.h"

namespace KetCat
{
	/// @brief Functor to generate a Gaussian wave packet state vector.
	/// @tparam Dim The dimension of the state vector to generate.
	/// @param x0     Center position
	/// @param k0     Central wave number 
	/// @param sigma  The standard deviation
	/// @param dx     Discretisation step
	template<dimension_t Dim>
	struct FreeParticleGaussianWavePacket
	{
		constexpr StateVector<Dim> operator()(float_t x0, float_t k0, float_t sigma, float_t dx) const noexcept
		{
			StateVector<Dim> psi = {};

			for (dimension_t n = 0; n < Dim; ++n)
			{
				// Position corresponding to index n
				const float_t x = (n + 1) * dx;

				// Gaussian envelope calculation 
				// exp(-((x - x0)^2) / (4 * sigma^2))
				const float_t exponent = -((x - x0) * (x - x0)) / (4.0 * sigma * sigma);
				const float_t envelope = ConstexprMath::exp<20>(exponent);

				// Plane wave component calculation: cos(k0 * x) + i * sin(k0 * x)
				const float_t realPart = ConstexprMath::cos(k0 * x);
				const float_t imagPart = ConstexprMath::sin(k0 * x);

				// Combine envelope and plane wave to form the complex amplitude
				psi[n] = cplx_t(envelope * realPart, envelope * imagPart);
			}

			//psi.normalize();
			return psi;
		}
	};
}