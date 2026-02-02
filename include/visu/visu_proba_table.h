#pragma once
#include <iostream>
#include "hilbert_space/state_vector.h"


namespace KetCat::Visu
{
   	/// @brief Simple terminal-based probability table visualization implementation
	/// @tparam Dim Dimension of the state vector
	/// @detail Displays probability densities as a table of numerical values.
	template<dimension_t Dim>
	struct VisuProbaTable
	{
		/// @brief Print the measurement probabilities for the selected qubits.
		/// @tparam SelectedQBits  Indices of the qubits to include in the probability table. 
        template<dimension_t... SelectedQBits>
        void update(const StateVector<FiniteHilbertSpace<Dim>> s) const
        {
            constexpr dimension_t NumSelected = sizeof...(SelectedQBits);
            constexpr qbit_list_t<NumSelected> SelectedQubits =
                qbit_list_t<NumSelected>{ static_cast<dimension_t>(SelectedQBits)... };

            // Format numbers
            std::cout << std::fixed << std::setprecision(2);

            // Number of qubits we are considering
            constexpr dimension_t ReducedDim = ConstexprMath::pow2(NumSelected);

   			// Get full state probabilities
            probability_vector_t<Dim> probabilities = s.getProbabilities();

            // Array to accumulate probabilities for each reduced basis state
            probability_vector_t<ReducedDim> ReducedProbabilities = {};

            // Sum probabilities over all other qubits
            for (dimension_t i = 0; i < Dim; ++i)
            {
                float_t p = probabilities[i];

                dimension_t ReducedIdx = 0;
                for (dimension_t b = 0; b < NumSelected; ++b)
                {
                    const dimension_t q = SelectedQubits[b];
                    if (i & (1ULL << q))
                        ReducedIdx |= (1ULL << b);
                }

                ReducedProbabilities[ReducedIdx] += p;
            }

            // Print header
            std::cout << "| Binary | Decimal | Probability (%) |\n";
            std::cout << "|--------|---------|----------------|\n";

            // Print reduced probabilities
            for (dimension_t i = 0; i < ReducedDim; ++i)
            {
                // Binary string
                std::cout << "|";
                for (dimension_t b = NumSelected; b-- > 0;)
                    std::cout << ((i & (1ULL << b)) ? '1' : '0');
                std::cout << "|";

                // Decimal
                std::cout << " " << i << " |";

                // Probability
                std::cout << " " << (ReducedProbabilities[i] * 100.0) << " % |\n";
            }
        }
	};
} // namespace KetCat::Visu