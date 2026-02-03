#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"

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
		constexpr StateVector<InfiniteHilbertSpace<Dim>>
			operator()(float_t x0, float_t k0, float_t sigma, float_t dx) const noexcept
		{
			StateVector<InfiniteHilbertSpace<Dim>> Psi = {};

			for (dimension_t n = 0; n < Dim; ++n)
			{
				// Position corresponding to index n
				const float_t x = (n + 1) * dx;

				// Gaussian envelope calculation 
				// exp(-((x - x₀)²) / (4 * σ²))
				const float_t Exponent = -((x - x0) * (x - x0)) / (4.0 * sigma * sigma);
				const float_t Envelope = ConstexprMath::exp<20>(Exponent);

				// Plane wave component calculation: cos(k₀ * x) + i * sin(k₀ * x)
				const float_t RealPart = ConstexprMath::cos(k0 * x);
				const float_t ImagPart = ConstexprMath::sin(k0 * x);

				// Combine envelope and plane wave to form the complex amplitude
				Psi[n] = cplx_t(Envelope * RealPart, Envelope * ImagPart);
			}

			Psi.normalize(dx);
			return Psi;
		}
	};
}