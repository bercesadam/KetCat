#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "quantum_number.h"


namespace KetCat
{
	/// @brief Helper function which computes the associated Laguerre polynomial L_p^(α)(x).
	/// @param p     Degree of the polynomial (non-negative integer).
	/// @param alpha Parameter of the polynomial (non-negative integer).
	/// @param x     Point at which to evaluate the polynomial.
	/// @return Value of the associated Laguerre polynomial L_p^(α)(x).
	static constexpr double laguerre(unsigned p, unsigned alpha, double x) noexcept
	{
		if (p == 0)
		{
			return 1.0;
		}

		double Lkm1 = 1.0;           // L_0^(α)(x)
		double Lk = 1.0 + alpha - x; // L_1^(α)(x)

		for (unsigned k = 1; k < p; ++k)
		{
			const double a = (2.0 * k + 1.0 + alpha - x);
			const double b = (k + alpha);
			const double Lk1 = (a * Lk - b * Lkm1) / (k + 1.0);
			Lkm1 = Lk;
			Lk = Lk1;
		}
		
		return Lk;
	}


	/// @brief  Construct a hydrogenic-like reduced radial wavefunction seed u(r) flattened to 1D.
	/// @param  n     Principal quantum number n ≥ 1.
	/// @param  l     Orbital angular momentum ℓ with 0 ≤ ℓ < n.
	/// @param  a_eff Effective length scale a_eff > 0.
	/// @param  dx    Grid spacing Δr.
	/// @return       Reduced radial component u(r), normalized so that Σ |u|² · Δr = 1.
	/// @details
	///   This returns u(r) (1D) which depends only on (n, ℓ). The full wavefunction is
	///   ψ_{nℓm}(r,θ,φ) = (u_{nℓ}(r)/r) · Y_{ℓm}(θ,φ). The Hamiltonian is m-independent
	///   for central potentials; m enters only via the angular factor Y_{ℓm}.
	///
	/// @tparam Dim Size of the discrete spatial grid
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
	struct HydrogenOrbital
	{
		/// @brief Generates a hydrogen-like orbital wavefunction.
		///
		/// @param q   Quantum number pair (n, l)
		/// @param a0  Effective Bohr radius (controls spatial scale)
		/// @param dx  Spatial discretization step
		/// @param x0  Position of the atomic center
		///
		/// @return Normalized quantum state vector representing the orbital
		template <quantum_number_t QuantumNumberType>
		constexpr StateVector<HilbertSpace>
			operator()(QuantumNumberType q, double a_eff) const noexcept
		{
			const natural_t n = q.n();
			const natural_t l = q.l();

			StateVector<HilbertSpace> Psi{ complex_t::zero() };

			// Radial grid: r_i = i·dx, i = 0..Dim−1; u(0) remains 0
			for (natural_t i = 1; i < HilbertSpace::Dim; ++i)
			{
				const double r = i * HilbertSpace::dx;
				const double x = 2.0 * r / (n * a_eff);

				// Behavior near r=0 (the nucleus)
				// r^(ℓ+1)
				double rPow = 1.0;
				for (unsigned k = 0; k < l + 1; ++k)
				{
					rPow *= r;
				}

				// Exponential tail
				// exp(−r / (n·a_eff))
				const double Exponential = ConstexprMath::exp<30>(-r / (n * a_eff));

				// Associated Laguerre: L_{n−ℓ−1}^(2ℓ+1)(x)
				const unsigned p = n - l - 1;
				const unsigned alpha = 2 * l + 1;
				const double Laguerre = laguerre(p, alpha, x);

				const double Value = rPow * Exponential * Laguerre;
				Psi[i] = complex_t::fromReal(Value);
			}

			// Enforce discrete radial normalization: Σ |u|² · Δr = 1
			Psi.normalize();

			return Psi;
		}
	};
}
