#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"
#include "atomic_physics_core/atom.h"
#include "atomic_physics_core/quantum_number.h"
#include "atomic_physics_core/rydberg_quantum_defect.h"


namespace KetCat
{
	/// @brief Helper function which computes the associated Laguerre polynomial L_p^(α)(x).
	///		   Currently unused legacy function, replaced by kummerHypergeometric1F1 to support non-integer principal quantum numebrs
	/// @param p     Degree of the polynomial (non-negative integer).
	/// @param alpha Parameter of the polynomial (non-negative integer).
	/// @param x     Point at which to evaluate the polynomial.
	/// @return Value of the associated Laguerre polynomial L_p^(α)(x).
	static constexpr real_t laguerre(unsigned p, unsigned alpha, real_t x) noexcept
	{
		if (p == 0)
		{
			return 1.0;
		}

		real_t Lkm1 = 1.0;           // L_0^(α)(x)
		real_t Lk = 1.0 + alpha - x; // L_1^(α)(x)

		for (unsigned k = 1; k < p; ++k)
		{
			const real_t a = (2.0 * k + 1.0 + alpha - x);
			const real_t b = (k + alpha);
			const real_t Lk1 = (a * Lk - b * Lkm1) / (k + 1.0);
			Lkm1 = Lk;
			Lk = Lk1;
		}
		
		return Lk;
	}

	/// @brief Computes the confluent hypergeometric function 1F1(a, b, x) (Kummer's function).
    /// @details Evaluates the series expansion 1 + (a/b)x + (a(a+1)/b(b+1))(x^2/2!) + ...
    /// This generalizes the associated Laguerre polynomials to non-integer degrees,
    /// which is required for calculating alkali Rydberg states with quantum defects.
	///
    /// @param a 	First parameter (numerator of the Pochhammer symbols).
    /// @param b 	Second parameter (denominator of the Pochhammer symbols).
    /// @param x 	The value at which to evaluate the function.
    /// @return The value of the Kummer confluent hypergeometric function.
	static constexpr real_t kummerHypergeometric1F1(real_t a, real_t b, real_t x) noexcept
	{
		real_t Term = 1.0;
		real_t Sum = 1.0;

		// Constrained to maximum 1000 terms as a safety limit
		for (natural_t k = 1; k < 1000; ++k)
		{
			Term *= (a + k - 1.0) / (b + k - 1.0) * (x / k);
			Sum += Term;

			// Stopping condition if the term is negligible
			if (ConstexprMath::abs(Term) < 1e-14 * ConstexprMath::abs(Sum)) 
				break;
		}
		return Sum;
	}
	
	/// @brief Construct a generalized hydrogenic reduced radial wavefunction seed u(r).
    ///
    /// This implements a radial bound-state orbital for a central Coulomb potential,
    /// generalized to support both pure hydrogenic states and alkali-metal Rydberg 
    /// states through the use of the quantum defect theory. 
    ///
    /// The returned function is the reduced radial wavefunction:
    ///
    ///     u(r) = r · R(r),
    ///
    /// discretized on a 1D radial grid and normalized on that grid.
    ///
    /// Unlike standard hydrogenic orbitals that rely on integer principal quantum 
    /// numbers (n), this implementation utilizes the effective principal quantum 
    /// number n* = n - δ_l, where δ_l is the l-dependent Rydberg quantum defect.
    ///
    /// To support non-integer values of n*, the associated Laguerre polynomial 
    /// L_{p}^{α}(x) is generalized via the Kummer confluent hypergeometric 
    /// function ₁F₁(a, b, x). This ensures that the nodal structure (radial nodes) 
    /// and the phase of the wavefunction correctly reflect the core penetration 
    /// effects in alkali atoms.
    ///
    /// The reduced radial function takes the form:
    ///
    ///     u_{n*l}(r) ∝ r^{l+1} · exp(−r / (n* a_eff)) · ₁F₁( -(n* - l - 1), 2l + 2, 2r / (n* a_eff) )
    ///
    /// When n* is an integer (e.g., pure Hydrogen), the Kummer series terminates 
    /// and the expression becomes mathematically identical to the standard 
    /// associated Laguerre polynomial form:
    ///
    ///     u_{nl}(r) ∝ r^{l+1} · exp(−r / (n a_eff)) · L_{n−l−1}^{2l+1}(2r / (n a_eff))
    ///
    /// The full spatial wavefunction is:
    ///
    ///     ψ_{n*lm}(r,θ,φ) = (u_{n*l}(r) / r) · Y_{lm}(θ,φ)
    ///
    /// @tparam HilbertSpace
    ///   Discrete 1D spatial Hilbert space defining the radial grid size and spacing.
    /// @tparam element
    ///   The chemical element used to retrieve the appropriate quantum defect and 
    ///   effective Bohr radius.
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace, Element element>
	struct HydrogenOrbitalRadial
	{
		/// @brief Generates a reduced radial part of a hydrogen-like orbital wavefunction.
		///
		/// @param element   Chemical element descriptor (e.g. Li, Na, K)
		/// @param q         Quantum number pair (n, l)
		///
		/// @return Normalized quantum state vector representing the radial part of a hydrogenic orbital
		template <quantum_number_t QuantumNumberType>
		constexpr StateVector<HilbertSpace>
			operator()(QuantumNumberType q) const noexcept
		{
			const real_t N_star = q.n() - RydbergQuantumDefect::value(element, q);
			const natural_t l = q.l();
			const real_t A_eff = Atom<element>::getEffectiveBohrRadius();

			StateVector<HilbertSpace> Psi{ complex_t::zero() };

			// Radial grid: r_i = i·dx, i = 0..Dim−1; u(0) remains 0
			for (natural_t i = 1; i < HilbertSpace::Dim; ++i)
			{
				const real_t r = HilbertSpace::gridToR(i);
				const real_t x = 2.0 * r / (N_star * A_eff);

				// Behavior near r=0 (the nucleus)
				// r^(ℓ+1)
				real_t rPow = 1.0;
				for (natural_t k = 0; k < l + 1; ++k)
				{
					rPow *= r;
				}

				// Exponential tail
				// exp(−r / (n·a_eff))
				const real_t Exponential = ConstexprMath::exp(-r / (N_star * A_eff));

				// Associated Laguerre: L_{n−ℓ−1}^(2ℓ+1)(x)
				const natural_t p = N_star - l - 1;
				const natural_t alpha = 2 * l + 1;
				const real_t Laguerre = laguerre(p, alpha, x);

				const real_t Value = rPow * Exponential * Laguerre;
				Psi[i] = complex_t::fromReal(Value);
			}

			// Enforce discrete radial normalization: Σ |u|² · Δr = 1
			Psi.normalize();

			return Psi;
		}
	};
}
