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
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
	struct FreeParticleGaussianWavePacket
	{
		constexpr StateVector<HilbertSpace>
			operator()(real_t x0, real_t k0, real_t sigma) const noexcept
		{
			constexpr real_t dx = HilbertSpace::dx;
			constexpr natural_t Dim = HilbertSpace::Dim;

			StateVector<HilbertSpace> Psi = {};

			for (natural_t n = 0; n < Dim; ++n)
			{
				// Position corresponding to index n
				const real_t x = n * dx;

				// Gaussian envelope calculation 
				// exp(-((x - x₀)²) / (4 * σ²))
				const real_t Exponent = -((x - x0) * (x - x0)) / (4.0 * sigma * sigma);
				const real_t Envelope = ConstexprMath::exp<20>(Exponent);

				// Plane wave component calculation: cos(k₀ * x) + i * sin(k₀ * x)
				const real_t RealPart = ConstexprMath::cos(k0 * x);
				const real_t ImagPart = ConstexprMath::sin(k0 * x);

				// Combine envelope and plane wave to form the complex amplitude
				Psi[n] = complex_t(Envelope * RealPart, Envelope * ImagPart);
			}

			// Dirichlet enforcement
			Psi[0] = complex_t::zero();
			Psi[Dim - 1] = complex_t::zero();

			Psi.normalize();
			return Psi;
		}
	};
}