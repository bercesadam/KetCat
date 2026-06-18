#pragma once
#include <string>
#include "atomic_units.h"
#include "laser/laser_pulse.h"
#include "hilbert_space/state_vector.h"


namespace KetCat
{
	/// @brief Data structure representing the output state of a single qubit (atom)
    /// for visualization and analysis purposes.
    template <hilbert_space_t FullHilbertSpace>
	struct QubitData
	{
		/// @brief Complex amplitudes of the qubit's state in the logical basis {|0⟩, |1⟩}
		complex_t m_alpha;
		complex_t m_beta;

		/// @brief Purity of the qubit's reduced density matrix (Tr(ρ²))
		real_t m_purity;

		/// @brief 2D spatial wavefunction of the atom corresponding to this qubit
		StateVector<FullHilbertSpace, QuantumPicture::Schrodinger> m_psi2D;

		/// @brief Laser parameters for this qubit (Wavelength and Intensity)
		real_t m_laser1Wavelength;
		real_t m_laser1Intensity;
		real_t m_laser2Wavelength;
		real_t m_laser2Intensity;

		/// @brief Custom caption for the simulation view
		std::string m_title;
	};


	/// @brief Data structure representing the complete state of the QPU at a given simulation time,
	/// including the global state vector and per-qubit information for visualization.
	template <hilbert_space_t FullHilbertSpace, natural_t QubitCount>
    struct SimulationView
    {
        /// @brief Simulation time in SI units  
        real_t m_time;

		/// @brief Global state vector, reduced to the logical subspace of the qubits
        probability_vector_t<ConstexprMath::pow(natural_t(2), QubitCount)> m_outputProbabilities;

		/// @brief Individual qubit data for each atom
        std::array<QubitData<FullHilbertSpace>, QubitCount> m_qubitDatum;

        /// @brief Custom caption for the simulation view -   TODO decide if needed 
        //std::string m_title;
    };
}
