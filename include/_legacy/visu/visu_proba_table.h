#pragma once
#include <iostream>
#include "hilbert_space/state_vector.h"


namespace KetCat::Visu
{
   	/// @brief Simple terminal-based probability table visualization implementation
	/// @tparam Dim Dimension of the state vector
	/// @detail Displays probability densities as a table of numerical values.
	template<natural_t Dim>
	struct VisuProbaTable
	{
		/// @brief Print the measurement probabilities for the selected qubits.
		/// @tparam SelectedQBits  Indices of the qubits to include in the probability table. 
        template<natural_t... SelectedQBits>
        void update(const StateVector<FiniteHilbertSpace<Dim>> s) const
        {
            constexpr natural_t NumSelected = sizeof...(SelectedQBits);
            constexpr qbit_list_t<NumSelected> SelectedQubits =
                qbit_list_t<NumSelected>{ static_cast<natural_t>(SelectedQBits)... };

            // Format numbers
            std::cout << std::fixed << std::setprecision(2);

            // Number of qubits we are considering
            constexpr natural_t ReducedDim = ConstexprMath::pow2(NumSelected);

   			// Get full state probabilities
            probability_vector_t<Dim> probabilities = s.getProbabilities();

            // Array to accumulate probabilities for each reduced basis state
            probability_vector_t<ReducedDim> ReducedProbabilities = {};

            // Sum probabilities over all other qubits
            for (natural_t i = 0; i < Dim; ++i)
            {
                real_t p = probabilities[i];

                natural_t ReducedIdx = 0;
                for (natural_t b = 0; b < NumSelected; ++b)
                {
                    const natural_t q = SelectedQubits[b];
                    if (i & (1ULL << q))
                        ReducedIdx |= (1ULL << b);
                }

                ReducedProbabilities[ReducedIdx] += p;
            }

            
            // ---- Column widths ----
            constexpr int BIN_WIDTH = NumSelected + 2; // padding inside cell
            constexpr int DEC_WIDTH = 8;
            constexpr int PROB_WIDTH = 14;

            auto line = []()
            {
                std::cout
                    << "+" << std::string(BIN_WIDTH, '-')
                    << "+" << std::string(DEC_WIDTH, '-')
                    << "+" << std::string(PROB_WIDTH, '-')
                    << "+\n";
            };

            // Print header
            line();
            std::cout << "| " << std::left << std::setw(BIN_WIDTH - 1) << "Bin"
                    << "| " << std::left << std::setw(DEC_WIDTH - 1) << "Dec"
                    << "| " << std::left << std::setw(PROB_WIDTH - 1) << "Proba (%)"
                    << "|\n";
            line();

            // Print rows
            for (natural_t i = 0; i < ReducedDim; ++i)
            {
                std::string bin;
                for (natural_t b = NumSelected; b-- > 0;)
                    bin += (i & (1ULL << b)) ? '1' : '0';

                std::cout << "| " << std::setw(BIN_WIDTH - 1) << std::right << bin
                        << "| " << std::setw(DEC_WIDTH - 1) << std::right << i
                        << "| " << std::setw(PROB_WIDTH - 1) << std::right
                        << std::fixed << std::setprecision(2)
                        << (ReducedProbabilities[i] * 100.0)
                        << " |\n";
            }

            line();

        }
	};
} // namespace KetCat::Visu