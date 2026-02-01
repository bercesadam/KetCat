#pragma once
#include "constexprmath/constexpr_core_functions.h"

namespace KetCat
{
   /// @brief One-dimensional Quantum Harmonic Oscillator potential.
   ///
   /// This class represents the standard 1D harmonic oscillator potential:
   ///
   /// V(x) = ½ m ω² (x - x₀)²
   ///
   /// where:
   /// - m : particle mass
   /// - ω : oscillator frequency
   /// - x₀: center of the oscillator (default 0)
	class HarmonicOscillatorPotential
	{
		/// Mass of the particle
		float_t m_m;
		/// Angular frequency ω
		float_t m_omega;
		/// Center position x₀
		float_t m_x0;

	public:
		/// @brief Constructs a harmonic oscillator potential.
		///
		/// @param m     Particle mass (default 1.0)
		/// @param omega Oscillator frequency ω (default 1.0)
		/// @param x0    Center position x₀ (default 0.0)
		constexpr HarmonicOscillatorPotential(float_t m,
			float_t omega, float_t x0) noexcept
			: m_m(m), m_omega(omega), m_x0(x0)
		{
		}

		/// @brief Evaluate the potential at a given position.
		///
		/// @param x Spatial coordinate
		/// @return V(x) = ½ m ω² (x - x₀)²
		constexpr float_t operator()(float_t x) const noexcept
		{
			const float_t dx = x - m_x0;
			return 0.5 * m_m * m_omega * m_omega * dx * dx;
		}
	};
}
