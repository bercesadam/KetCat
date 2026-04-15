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

	
	/// @brief Construct a hydrogenic reduced radial wavefunction seed u(r)
	///   This implements a hydrogen-like bound-state orbital for a Coulomb potential.
	///   The returned function is the reduced radial wavefunction
	///
	///     u(r) = r · R(r),
	///
	///   discretized on a 1D radial grid and normalized on that grid.
	///
	///   Unlike Slater-type orbitals, hydrogenic orbitals are exact eigenfunctions of the
	///   hydrogen atom Hamiltonian and are expressed in terms of associated Laguerre
	///   polynomials. The radial dependence is controlled by the principal quantum number n
	///   and orbital angular momentum ℓ.
	///
	///   The reduced radial function takes the form:
	///
	///     uₙℓ(r) ∝ r^{ℓ+1} · exp(−r / (n a_eff)) · L_{n−ℓ−1}^{2ℓ+1}(2r / (n a_eff))
	///
	///   where L_{p}^{α}(x) is an associated Laguerre polynomial and a_eff is an effective
	///   Bohr radius (e.g. incorporating screening or reduced-mass effects).
	///
	///   The full spatial wavefunction is
	///
	///     ψ_{nℓm}(r,θ,φ) = (u_{nℓ}(r) / r) · Y_{ℓm}(θ,φ)
	///
	///   As with all central potentials, the Hamiltonian is m-independent; the magnetic
	///   quantum number m enters only through the spherical harmonic Y_{ℓm}.
	///
	/// @tparam HilbertSpace
	///   Discrete 1D spatial Hilbert space defining the radial grid size and spacing.

	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
	struct HydrogenOrbital
	{
		/// @brief Generates a reduced radial part of a hydrogen-like orbital wavefunction.
		///
		/// @param element   Chemical element descriptor (e.g. Li, Na, K)
		/// @param q         Quantum number pair (n, l)
		///
		/// @return Normalized quantum state vector representing the radial part of a hydrogenic orbital
		template <quantum_number_t QuantumNumberType>
		constexpr StateVector<HilbertSpace>
			operator()(Element element, QuantumNumberType q) const noexcept
		{
			const natural_t n = q.n();
			const natural_t l = q.l();
			const real_t A_eff = Atom<element>::getEffectiveBohrRadius();

			StateVector<HilbertSpace> Psi{ complex_t::zero() };

			// Radial grid: r_i = i·dx, i = 0..Dim−1; u(0) remains 0
			for (natural_t i = 1; i < HilbertSpace::Dim; ++i)
			{
				const double r = i * HilbertSpace::dx;
				const double x = 2.0 * r / (n * A_eff);

				// Behavior near r=0 (the nucleus)
				// r^(ℓ+1)
				double rPow = 1.0;
				for (unsigned k = 0; k < l + 1; ++k)
				{
					rPow *= r;
				}

				// Exponential tail
				// exp(−r / (n·a_eff))
				const double Exponential = ConstexprMath::exp<30>(-r / (n * A_eff));

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
