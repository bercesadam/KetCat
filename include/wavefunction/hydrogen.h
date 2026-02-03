#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"

namespace KetCat
{
	/// @brief Represents a pair of quantum numbers (n, l) using the spectroscopic notation:
	/// - n : Principal quantum number (n >= 1), controls energy and radial extent
	/// - l : Orbital angular momentum quantum number (0 <= l <= n-1),
	///       controls angular structure and parity
	///	This is introduced to disallow using non-existing configurations in 
	/// atomic simulations and also as a convenience layer.
	class QuantumNumber
	{
		/// Principal quantum number
		unsigned int m_n;
		/// Orbital angular momentum quantum number
		unsigned int m_l;

		/// @brief Private constructor enforcing valid factory usage.
		/// @param n Principal quantum number
		/// @param l Orbital angular momentum quantum number
		constexpr QuantumNumber(unsigned int n, unsigned int l)
			: m_n(n), m_l(l)
		{
		}

	public:
		/// @brief Returns the principal quantum number.
		constexpr unsigned int n() { return m_n; }

		/// @brief Returns the orbital angular momentum quantum number.
		constexpr unsigned int l() { return m_l; }

		/// @name Hydrogen orbital presets
		/// @{
		static constexpr QuantumNumber _1s() { return QuantumNumber{ 1, 0 }; }
		static constexpr QuantumNumber _2s() { return QuantumNumber{ 2, 0 }; }
		static constexpr QuantumNumber _2p() { return QuantumNumber{ 2, 1 }; }
		static constexpr QuantumNumber _3s() { return QuantumNumber{ 3, 0 }; }
		static constexpr QuantumNumber _3p() { return QuantumNumber{ 3, 1 }; }
		static constexpr QuantumNumber _3d() { return QuantumNumber{ 3, 2 }; }
		static constexpr QuantumNumber _4s() { return QuantumNumber{ 4, 0 }; }
		static constexpr QuantumNumber _4p() { return QuantumNumber{ 4, 1 }; }
		static constexpr QuantumNumber _4d() { return QuantumNumber{ 4, 2 }; }
		static constexpr QuantumNumber _4f() { return QuantumNumber{ 4, 3 }; }
		/// @}
	};


	/// @brief Helper function which computes the associated Laguerre polynomial L_p^(α)(x).
	/// @param p     Degree of the polynomial (non-negative integer).
	/// @param alpha Parameter of the polynomial (non-negative integer).
	/// @param x     Point at which to evaluate the polynomial.
	/// @return Value of the associated Laguerre polynomial L_p^(α)(x).
	static constexpr double laguerre(unsigned p, unsigned alpha, double x) noexcept
	{
		if (p == 0) return 1.0;

		double Lkm1 = 1.0;                   // L_0^(α)(x)
		double Lk = 1.0 + alpha - x;       // L_1^(α)(x)

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
	template<dimension_t Dim>
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
		constexpr StateVector<InfiniteHilbertSpace<Dim>>
			operator()(QuantumNumber q, double a_eff, double dx) const noexcept
		{
			const unsigned int n = q.n();
			const unsigned int l = q.l();

			StateVector<InfiniteHilbertSpace<Dim>> Psi{ cplx_t::zero() };

			// Radial grid: r_i = i·dx, i = 0..Dim−1; u(0) remains 0
			for (dimension_t i = 1; i < Dim; ++i)
			{
				const double r = i * dx;
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
				Psi[i] = cplx_t::fromReal(Value);
			}

			// Enforce discrete radial normalization: Σ |u|² · Δr = 1
			Psi.normalize(dx);

			return Psi;
		}
	};
}
