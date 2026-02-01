#pragma once
#include "constexprmath/constexpr_core_functions.h"

namespace KetCat
{
	/// @brief Radial soft-Coulomb potential functor 
	/// @details
	/// This class implements a radial soft-Coulomb potential commonly used
	/// in  simulations of atomic systems. The potential includes both
	/// the Coulomb attraction term and the centrifugal barrier term for
	/// angular momentum states.
	class SoftCoulombRadialPotential
	{
		double m_Zeff;     // Effective nuclear charge Z_eff
		double m_a;        // Softening parameter a > 0
		unsigned m_l;      // Orbital quantum number ℓ ≥ 0
		double m_hbar;     // ℏ
		double m_mu;       // Reduced mass μ

	public:
		/// @brief Constructs a radial soft-Coulomb potential.
		/// @param  Zeff   Effective nuclear charge Z_eff
		/// @param  a      Softening parameter a > 0
		/// @param  l      Orbital quantum number ℓ ≥ 0
		/// @param  hbar   Reduced Planck's constant ℏ
		/// @param  mu     Reduced mass μ
		constexpr SoftCoulombRadialPotential(double Zeff,
			double a, unsigned l, double hbar, double mu) noexcept
			: m_Zeff(Zeff), m_a(a), m_l(l), m_hbar(hbar), m_mu(mu)
		{
		}

		/// @brief  Evaluate V(r) at radius r.
		/// @param  r   Position from the nuclear center.
		/// @return     Potential value V(r).
		///
		/// @details
		/// Computation:
		///   r² + a² → r2s
		///   V(r) = −Z_eff / √(r2s)  +  ℓ(ℓ+1)·ℏ² / (2μ·r2s)
		constexpr double operator()(double r) const noexcept
		{
			const double SquareRadiusWithSoftening = r * r + m_a * m_a;  // r² + a²
			const double CoulombAttraction = -m_Zeff / ConstexprMath::sqrt(SquareRadiusWithSoftening);
			const double CentrifugalBarrier = (static_cast<double>(m_l) * (m_l + 1) * m_hbar * m_hbar) / (2.0 * m_mu * SquareRadiusWithSoftening);
			return CoulombAttraction + CentrifugalBarrier;
		}
	};
}
