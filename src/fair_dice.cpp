#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace KetCat::QCC;

///@brief Fair dice algorithm using 3 qubits
///       Demonstrating how rotations on the Bloch sphere can be eningeered
///       to achieve a desired probability distribution and also demonstrates
///       that quantum gates are always reversible (uncomputing)

int main()
{
	std::cout << "Fair dice from 3 qubits\n";

	// Choose θ such that after Ry(θ):
	//   P(|1>) = sin²(θ/2) = 1/3
	//   P(|0>) = 2/3
	//   θ = 2 * arcsin(√(1/3)) ≈ 1.23096
	constexpr auto Ry13_23 = Gates::RotationY(1.23095941734);

	constexpr auto RyPi4 = Gates::RotationY(ConstexprMath::Pi / 4);
	constexpr auto RyMinusPi4 = Gates::RotationY(-ConstexprMath::Pi / 4);

	constexpr auto FairDice = QuantumCircuit<3>().withGates(
		// Prepare qubit 2 with biased probabilities:
		//   P(|1>) = 1/3, P(|0>) = 2/3
		QuantumGate<1, Ry13_23>().toBits(2),

		// Flip qubit 2 so that:
		//   P(|1>) = 2/3, P(|0>) = 1/3
		// This makes |1> the "common" branch
		QuantumGate<1, Gates::X>().toBits(2),

		// Rotate qubit 1 into a |0>/|1> superposition
		// Used as a workspace qubit to split probability mass
		QuantumGate<1, RyPi4>().toBits(1),

		// Entangle qubit 1 with qubit 2
		// This redistributes probability conditioned on qubit 1
		QuantumGate<2, Gates::CX>().toBits(2, 1),

		// Undo the rotation on qubit 1 (uncompute workspace state)
		QuantumGate<1, RyMinusPi4>().toBits(1),

		// Remove entanglement between qubits 2 and 1
		// In this current program this has no real effects
		// but in real quantum computers and simulators
		// this would be crucial to make the measurements the these
	 	// two bits independent again
		QuantumGate<2, Gates::CX>().toBits(2, 1),

		// Flip qubit 1 to align the probability structure
		// so that each (qubit2, qubit1) pair has total weight 1/3
		QuantumGate<1, Gates::X>().toBits(1),

		// Apply Hadamard on qubit 0 to split each 1/3 branch into two
		// → 3 × 2 = 6 equiprobable outcomes (fair dice)
		QuantumGate<1, Gates::H>().toBits(0)
	);

    KetCat::Visu::VisuProbaTable<8>().update<0, 1, 2>(FairDice.getStateVector());

	return 0;
}

