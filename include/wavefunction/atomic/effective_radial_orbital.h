#pragma once
#include "hartree.h"
#include "hydrogenic_radial.h"
#include "slater_type_radial.h"

#include "atomic_physics_core/atom.h"
#include "atomic_physics_core/quantum_number.h"
#include "atomic_physics_core/rydberg_quantum_defect.h"

	
namespace KetCat
{
	/// @brief  This is a meta-wavefunction generator that performs a model selection for the radial part
    ///         for alkali atom 2D/3D wavefunction modeling
	///
	/// @details
	///   For alkali atoms, the valence electron experiences a Coulomb-like potential at long range, but the inner electrons cause significant screening and deviations from hydrogenic behavior.
	///   This generator uses quantum numbers and quantum defects to decide whether to use a hydrogenic orbital seed or a Slater-type orbital (STO) seed for the radial part of the wavefunction.
	///	  The decision is based on:
	///   - For high angular momentum (l >= 4) and sufficiently high principal quantum number (n >= 6), the valence electron is mostly outside the core and can be approximated by a hydrogenic orbital.
	///   - For d orbitals (l = 2) with n >= 6, if the quantum defect is small (< 0.1), it indicates near-hydrogenic behavior, so a hydrogenic seed is used.
	///   For all other cases, a Slater-type orbital seed is generated using an effective nuclear charge calculated from Slater's rules.
	///
	/// @tparam HilbertSpace 1D Spatial Hilbert space defining the radial grid
	template<spatial_hilbert_space_with_dim_t<1_D> HilbertSpace, Element element>
	struct EffectiveRadialOrbital
	{
		/// @brief Generates the effective radial orbital seed for alkali atoms,
        /// using either hydrogenic or Slater-type orbitals based on quantum numbers and quantum defects.
		///
		/// @param element   Chemical element descriptor (e.g. Li, Na, K)
		/// @param q         Quantum number pair
		///
		/// @return Normalized quantum state vector representing the orbital
		template <quantum_number_t QuantumNumberType>
		constexpr Wavefunction<HilbertSpace> operator()(QuantumNumberType q) const noexcept
		{
			const natural_t l = q.l();
			const real_t QuantumDefect = RydbergQuantumDefect::value(element, q);

			StateVector<HilbertSpace> Psi{ complex_t::zero() };
			
			bool UseHydrogenic = false;

			// Decision logic for model selection based on the azimuthal quantum number:
			// for f orbitals and higher (l >= 3), the valence electron is mostly outside the core and can be approximated by a hydrogenic orbital.
            if (l >= 3)
			{
				UseHydrogenic = true;
			}

			// In case of d orbitals, let's check the quantum defect value: if it's small, it indicates near-hydrogenic behavior, so we can use a hydrogenic seed.
			else if (l == 2 && QuantumDefect < 0.1)
			{
				UseHydrogenic = true;
			}

			// In case the quantum defect is negligible, we can also use a hydrogenic seed, even for lower l values.
			// This also effectively chooses the hydrogenic seed for actual Hydrogen, where all quantum defects are zero.
			else if (QuantumDefect < 0.05)
			{
				UseHydrogenic = true;
			}

			if (UseHydrogenic)
			{
				Psi = HydrogenOrbitalRadial<HilbertSpace, element>{}(q);
			}
			else
			{
				Psi = SlaterOrbitalRadial<HilbertSpace, element>{}(q);
			}

			constexpr real_t HartreeEnergy = calculateHartreeEnergy(element, q);
			return { Psi, HartreeEnergy };
		}
	};
}
