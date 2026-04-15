#pragma once
#include <utility>
#include "constants.h"
#include "elements.h"


namespace KetCat
{
	// Maximum number of shells/subshells to consider in the electron configuration
	constexpr natural_t MAX_SHELLS = 10; 

	// Type alias for the electron configuration array, which holds the electron counts for each shell/subshell.
	using electron_config_t = std::array<ElectronShell, MAX_SHELLS>;

	/// @brief Helper struct to store the number of electrons on each subshell
	struct ElectronShell
	{
		natural_t m_n; // Principal quantum number
		natural_t m_l; // Orbital angular momentum quantum number
		natural_t m_numElectrons; // Number of electrons in this shell
	};

	/// @brief Atom class template to store basic atomic information and electron configuration for seed wavefunction generation.
	/// @tparam E Element type (e.g. Element::Li, Element::Na)
	template <Element E>
	class Atom
	{
		// @brief Internal struct to hold the electron configuration and outer shell index for the atom.
		struct AtomData
		{
			// Atomic number (total number of electrons in a neutral atom)
			natural_t m_Z;
			
			// Effective Bohr radius
			real_t m_Aeff

			// Electron configuration array, where each entry corresponds to a subshell defined by (n, l)
			//and the number of electrons in that subshell.
			electron_config_t m_ElectronConfiguration;

			// Index of the outermost occupied shell/subshell in the electron configuration.
			// (Index of the last used element in the static array.)
			natural_t m_OuterShellIndex;
		};

		// @brief Electron configuration and outer shell index for the atom, generated at compile time.
		static constexpr AtomData m_Data = generateConfig();

		///@brief Calculate the electron configuration for a given element based on its atomic number Z.
		static constexpr AtomData generateConfig()
		{
			AtomData Data{};

			// Store the atomic number in the data struct
			Data.m_Z = AtomicNumber<E>::value;

			// Calculate effective Bohr radius
			// WARNING! Defining it as equal to the Bohr radius as currently this is used only
			// for the calculation of Hydrogenic like radial orbitals in Rydberg states,
			// where the effective Bohr radius is close to the actual Bohr radius.
			Data.m_Aeff = BohrRadius;

			// Temporary struct to represent a subshell with its quantum numbers (n, l).
			struct SubshellStub { natural_t n, l; };

			// The Aufbau principle dictates the order in which electrons fill atomic orbitals.
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
			natural_t Remaining = Z;

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

		/// @brief Get the effective Bohr radius for this atom, which can be used in hydrogenic orbital calculations.
		/// @return The effective Bohr radius (currently set equal to the actual Bohr radius for simplicity).
		static constexpr real_t getEffectiveBohrRadius() noexcept
		{
			return m_Aeff;
		}

		/// @brief Get the atomic number (total number of electrons) for this atom.
		/// @return The atomic number Z.
		static constexpr natural_t getAtomicNumber() noexcept
		{
			return m_Z;
		}
	};	
}
