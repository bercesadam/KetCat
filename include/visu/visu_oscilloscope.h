#pragma once
#include <chrono>
#include <thread>
#include <iostream>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <optional>
#include <functional>

#include "oscilloscope_config.h"
#include "hilbert_space/state_vector.h"
#include "constexprmath/constexpr_trigon.h"
#include "hamiltonian/hamiltonian.h" // for potential_functor concept
#include "hamiltonian/potential_barrier.h" // for ZeroPotential

#ifdef _WIN32
#include <windows.h>
//Funny easter egg in windows.h, so this is needed otherwise it breaks std::max
#undef max
#endif

namespace KetCat::Visu
{
	/// @brief Map phase angle arg(ψ) to ANSI color code
	inline const char* phaseToColor(real_t phase)
	{
		if (phase < -ConstexprMath::Pi * 0.5) return "\x1B[34m";  // dark blue
		if (phase < -ConstexprMath::Pi * 0.25) return "\x1B[94m"; // light blue
		if (phase < ConstexprMath::Pi * 0.25) return "\x1B[97m";  // white
		if (phase < ConstexprMath::Pi * 0.5) return "\x1B[91m";   // light red
		return "\x1B[31m";                                        // red
	}

	/// @brief Render a single oscilloscope line from (value, color) samples
	///
	/// @tparam Dim Dimension of the signal
	/// @param samples Array of tuples: (value, ANSI color)
	/// @param label Text label printed before the line
	template<dimension_t Dim>
	inline void renderLine(
		const std::array<std::tuple<real_t, const char*>, Dim>& samples,
		const char* label)
	{
		static const char* bars[] = {
			"\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83",
			"\xE2\x96\x84", "\xE2\x96\x85", "\xE2\x96\x86",
			"\xE2\x96\x87", "\xE2\x96\x88"
		};

		// Find maximum absolute value for normalization
		real_t maxVal = 0.0f;
		for (const auto& s : samples)
		{
			maxVal = std::max(maxVal, std::abs(std::get<0>(s)));
		}
		if (maxVal == 0.0f)
		{
			maxVal = 1.0f;
		}

		std::cout << label << " |";
		for (dimension_t i = 0; i < Dim; ++i)
		{
			const real_t value = std::get<0>(samples[i]);
			const char* color = std::get<1>(samples[i]);

			const real_t norm = std::abs(value) / maxVal;
			const std::size_t idx = static_cast<std::size_t>(norm * 7);

			std::cout << color << bars[idx] << "\x1B[0m";
		}
		std::cout << "|\n";
	}

	/// @brief Evaluate a potential functor over discrete spatial points
	/// @tparam Dim Dimension of the spatial discretization
	/// @param potential The potential functor to evaluate
	/// @param dx Spatial discretization step
	/// @return Array of potential values at discrete points
	template <dimension_t Dim, typename PotentialFunctor>
		requires potential_functor<PotentialFunctor, real_t>
	std::array<real_t, Dim> evalutePotentialFunctor(const PotentialFunctor& potential, const real_t dx)
	{
		std::array<real_t, Dim> DiscretePotentials{};
		for (dimension_t i = 0; i < Dim; ++i)
		{
			// Calculate position for the Potential callable: i * Δx
			const real_t x = i * dx;
			// Evaluate potential at position x
			DiscretePotentials[i] = potential(x);
		}
		return DiscretePotentials;
	}

	/// @brief Terminal-based oscilloscope visualization for 1D quantum states
	///
	/// @tparam Dim Dimension of the state vector
	///
	/// @details
	/// Displays:
	///  - Optional real part Re(ψ)  (yellow)
	///  - Optional imaginary part Im(ψ) (cyan)
	///  - Probability density |ψ|² (optionally phase-colored)
	template<dimension_t Dim>
	class VisuOscilloscope
	{
		/// Configuration options
		UsePhaseEncoding m_usePhaseEncoding;
		ClearScreen m_clearScreen;
		ShowComplexParts m_showComplex;
		ShowPotential m_showPotential;

		/// Optional parameters for potential visualization
		std::optional<real_t> m_dx;
		std::optional<std::function<real_t(real_t)>> m_potential;
		
	public:
		/// @brief Set the potential functor and spatial step for visualization
		/// @tparam PotentialFunctor Type of the potential functor
		/// @param potential The potential functor to visualize
		/// @param dx Spatial discretization step
		template <typename PotentialFunctor>
				requires potential_functor<PotentialFunctor, real_t>
		void setPotential(const PotentialFunctor& potential, real_t dx)
		{
			m_potential = potential;
			m_dx = dx;
		}

		/// @brief Clear the stored potential functor and spatial step
		void clearPotential()
		{
			m_potential.reset();
			m_dx.reset();
		}

		/// @brief Construct a VisuOscilloscope with specified settings
		/// @param usePhaseEncoding Enable phase-based coloring of probability density
		/// @param clearScreen      Clear screen before rendering
		/// @param showComplex      Enable visualization of real and imaginary parts
		/// @param showPotential    Enable visualization of potential V(x)
		VisuOscilloscope(
			const UsePhaseEncoding usePhaseEncoding,
			const ClearScreen clearScreen,
			const ShowComplexParts showComplex,
			const ShowPotential showPotential)
			: m_usePhaseEncoding(usePhaseEncoding),
			  m_clearScreen(clearScreen),
			  m_showComplex(showComplex),
			  m_showPotential(showPotential)
		{
#ifdef _WIN32
			// Enable UTF-8 output on Windows console
			SetConsoleOutputCP(CP_UTF8);
#endif
		}
		
		/// @brief Update the visualization with the current state vector
		/// @param s Current state vector
		template<KetCat::hilbert_space_t HilbertSpace>
		void update(const StateVector<HilbertSpace>& s) const
		{
			using namespace std::chrono_literals;

			if (enabled(m_clearScreen))
			{
				std::cout << "\x1B[2J\x1B[H";
			}

			// --- Probability density |ψ|² ---
			std::array<std::tuple<real_t, const char*>, Dim> ProbLine{};

			for (dimension_t i = 0; i < Dim; ++i)
			{
				const real_t p = s[i].normSquared();

				const char* color = "\x1B[97m";
				if (enabled(m_usePhaseEncoding))
				{
					const real_t Phase = std::atan2(s[i].im, s[i].re);
					color = phaseToColor(Phase);
				}

				ProbLine[i] = { p, color };
			}

			renderLine<Dim>(ProbLine, "Proba:     ");


			// --- Optional: Real and Imaginary parts ---
			if (enabled(m_showComplex))
			{
				std::array<std::tuple<real_t, const char*>, Dim> ReLine{};
				std::array<std::tuple<real_t, const char*>, Dim> ImLine{};

				for (dimension_t i = 0; i < Dim; ++i)
				{
					ReLine[i] = { s[i].re, "\x1B[33m" }; // Yellow
					ImLine[i] = { s[i].im, "\x1B[36m" }; // Cyan
				}

				renderLine<Dim>(ReLine, "Real:      ");
				renderLine<Dim>(ImLine, "Imag:      ");
			}

			// --- Optional: Potential V(x) ---
			if (enabled(m_showPotential))
			{
				if (!m_potential || !m_dx)
				{
                	std::cout<< "Warning! Ocilloscope is configured with ShowPotential=YES but potential/dx were never provided!" << std::endl;
                }

				const auto Potentials = evalutePotentialFunctor<Dim>(*m_potential, *m_dx);
				std::array<std::tuple<real_t, const char*>, Dim> PotentialLine{};
				for (dimension_t i = 0; i < Dim; ++i)
				{
					// Evaluate potential at position x
					PotentialLine[i] = { Potentials[i], "\x1B[35m"}; // Magenta
				}

				renderLine<Dim>(PotentialLine, "Potential: ");
			}

			// Small delay to allow visualization update
			std::this_thread::sleep_for(100ms);
		}
	};
}
