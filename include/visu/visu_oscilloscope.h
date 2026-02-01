#pragma once
#include <chrono>
#include <thread>
#include <iostream>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>

#include "wavefunction/state_vector.h"
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
	/// @brief Enum to specify whether to use phase encoding in visualization
	enum class UsePhaseEncoding
	{
		YES,
		NO
	};

	/// @brief Enum to specify whether to clear the screen before updating visualization
	enum class ClearScreen
	{
		YES,
		NO
	};

	/// @brief Enum to enable visualization of real and imaginary parts
	enum class ShowComplexParts
	{
		YES,
		NO
	};

	/// @brief Enum to specify whether to show potential in visualization
	enum class ShowPotential
	{
		YES,
		NO
	};

	/// @brief Check if an enum flag is enabled
	template<typename EnumType>
	constexpr inline bool enabled(EnumType e, EnumType yes = EnumType::YES) noexcept
	{
		return e == yes;
	}

	/// @brief Map phase angle arg(ψ) to ANSI color code
	inline const char* phaseToColor(float_t phase)
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
		const std::array<std::tuple<float_t, const char*>, Dim>& samples,
		const char* label)
	{
		static const char* bars[] = {
			"\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83",
			"\xE2\x96\x84", "\xE2\x96\x85", "\xE2\x96\x86",
			"\xE2\x96\x87", "\xE2\x96\x88"
		};

		// Find maximum absolute value for normalization
		float_t maxVal = 0.0f;
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
			const float_t value = std::get<0>(samples[i]);
			const char* color = std::get<1>(samples[i]);

			const float_t norm = std::abs(value) / maxVal;
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
		requires potential_functor<PotentialFunctor, float_t>
	std::array<float_t, Dim> evalutePotentialFunctor(const PotentialFunctor& potential, const float_t dx)
	{
		std::array<float_t, Dim> DiscretePotentials{};
		for (dimension_t i = 0; i < Dim; ++i)
		{
			// Calculate position for the Potential callable: i * Δx
			const float_t x = i * dx;
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
	struct VisuOscilloscope
	{
		/// Configuration options
		UsePhaseEncoding m_usePhaseEncoding;
		ClearScreen m_clearScreen;
		ShowComplexParts m_showComplex;
		ShowPotential m_showPotential;

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
		/// @param s Current quantum state vector
		/// @param usePhaseEncoding Enable phase-based coloring of probability density
		/// @param cls Clear screen before rendering
		/// @param showComplex Enable visualization of real and imaginary parts
		template<typename PotentialFunctor>
			requires potential_functor<PotentialFunctor, float_t>
		void update(const StateVector<Dim>& s, const PotentialFunctor potential, const float_t dx) const
		{
			using namespace std::chrono_literals;

			if (enabled(m_clearScreen))
			{
				std::cout << "\x1B[2J\x1B[H";
			}

			// --- Probability density |ψ|² ---
			std::array<std::tuple<float_t, const char*>, Dim> ProbLine{};

			for (dimension_t i = 0; i < Dim; ++i)
			{
				const float_t p = s[i].normSquared();

				const char* color = "\x1B[97m";
				if (enabled(m_usePhaseEncoding))
				{
					const float_t Phase = std::atan2(s[i].im, s[i].re);
					color = phaseToColor(Phase);
				}

				ProbLine[i] = { p, color };
			}

			renderLine<Dim>(ProbLine, "Proba:     ");


			// --- Optional: Real and Imaginary parts ---
			if (enabled(m_showComplex))
			{
				std::array<std::tuple<float_t, const char*>, Dim> ReLine{};
				std::array<std::tuple<float_t, const char*>, Dim> ImLine{};

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
				const auto Potentials = evalutePotentialFunctor<Dim>(potential, dx);
				std::array<std::tuple<float_t, const char*>, Dim> PotentialLine{};
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

		/// @brief Update the visualization with the current state vector (no potential)
		/// @param s Current quantum state vector
		/// @note This overload uses a zero potential by default
		void update(const StateVector<Dim>& s) const
		{
			update(s, ZeroPotential, 0.0);
		}
	};
}
