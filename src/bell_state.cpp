#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace KetCat::QCC;

//@brief Sanity check: Create and run a Bell Ψ⁺ state circuit on 2 qubits
//       Expected state: (|00⟩ + |11⟩) / √2

int main()
{
	std::cout << "Bell State Circuit demo (2 qubits)\n";

    constexpr auto BellStateCircuit = QuantumCircuit<2>().withGates(
        QuantumGate<1, Gates::H>().toBits(0),
        QuantumGate<2, Gates::CX>().toBits(0, 1));

    KetCat::Visu::VisuProbaTable<4>().update<0, 1>(BellStateCircuit.getStateVector());

    return 0;
}

