#pragma once
#include <utility>
#include "elements.h"

namespace KetCat
{
	// Maximum number of shells/subshells to consider in the electron configuration
	constexpr natural_t MAX_SHELLS = 10; 

	// Type alias for the electron configuration array, which holds the electron counts for each shell/subshell.
	using electron_config_t = std::array<ElectronShell, MAX_SHELLS>;

	/// @brief Helper struct to calculate 
	struct ElectronShell
	{
		natural_t m_n; // Principal quantum number
		natural_t m_l; // Orbital angular momentum quantum number
		natural_t m_numElectrons; // Number of electrons in this shell
	};

	/// @brief Atom class template that currently generates the electron configuration for a given element based on its atomic number.
	/// (Preserved for potential future use of more detailed atomic properties.)
	template <Element E>
	class Atom
	{
		// @brief Internal struct to hold the electron configuration and outer shell index for the atom.
		struct AtomData
		{
			electron_config_t m_ElectronConfiguration;
			natural_t m_OuterShellIndex;
		};

		// @brief Electron configuration and outer shell index for the atom, generated at compile time.
		static constexpr AtomData m_Data = generateConfig();

		///@brief Calculate the electron configuration for a given element based on its atomic number Z.
		static constexpr AtomData generateConfig()
		{
			AtomData Data{};

			struct SubshellStub { natural_t n, l; };
			constexpr SubshellStub Aufbau[] =
			{
				{1,0},
				{2,0}, {2,1},
				{3,0}, {3,1},
				{4,0}, {3,2}, {4,1},
				{5,0}, {4,2}, {5,1},
				{6,0}, {4,3}, {5,2},
				{6,1}, {7,0}, {5,3},
				{6,2}, {7,1},
				{8,0}, {5,4}, {6,3}, {7,2}, {8,1},
				{9,0}, {6,4}, {7,3}, {8,2}, {9,1},
				{10,0}
			};

			// The atomic number Z corresponds to the total number of electrons in a neutral atom,
			// which is equal to the underlying value of the Element enum.
			natural_t Remaining = std::to_underlying(E);

			for (natural_t i = 0; i < std::size(Aufbau); ++i)
			{
				if (Remaining <= 0)
				{
					break;
				}

				// Each subshell can hold a maximum of 2(2l + 1) electrons
				natural_t SubshellCapacity = 2 * (2 * Aufbau[i].l + 1);

				// Add electrons to the current subshell, but don't exceed the remaining electrons or the subshell capacity
				natural_t NoElectrons = (Remaining < SubshellCapacity) ? Remaining : SubshellCapacity;

				Data.m_ElectronConfiguration[i] = { Aufbau[i].n, Aufbau[i].l, NoElectrons };
				Remaining -= NoElectrons;
				Data.m_OuterShellIndex = i; 
			}

			return Data;
		}

	public:
		/// @brief Get the electron configuration for this atom.
		/// @return A constexpr array of ElectronShell structs representing the electron configuration.
		static constexpr electron_config_t getElectronConfiguration() noexcept
		{
			return m_ElectronConfiguration;
		}

		/// @brief Get the index of the outermost occupied shell/subshell in the electron configuration.
		/// @return The index of the outermost occupied shell/subshell.
		static constexpr natural_t getOuterShellIndex() noexcept
		{
			return m_OuterShellIndex;
		}
	};
}
