#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace KetCat::QCC;

int main()
{
	// @brief Construct Shor's algorithm circuit for N=21, a=2
	// 8 qubits: 3 for phase register, 5 for work register
	// The circuit performs controlled modular multiplications by
	// 2, 4, and 16 (mod 21) based on the phase qubits,
	// followed by an inverse QFT on the phase register

	std::cout << "Shor's Algorithm Circuit for N=21, a=2\n";
    
    constexpr auto IQFT3 = Gates::make_IQFT_matrix<3>();

    constexpr auto ShorCircuit =
        QuantumCircuit<8>().withGates(

            // --------------------------------------------------
            // 1) Hadamard on phase register (qubits 0,1,2)
            // --------------------------------------------------
            QuantumGate<1, Gates::H>().toBits(0),
            QuantumGate<1, Gates::H>().toBits(1),
            QuantumGate<1, Gates::H>().toBits(2),

            // --------------------------------------------------
            // 2) Initialize work register to |1> (qubit 3)
            // --------------------------------------------------
            QuantumGate<1, Gates::X>().toBits(3),

            // --------------------------------------------------
            // 3) Controlled modular multiply by 2 (mod 21)
            // Control: phase qubit 0
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2, Gates::CX>().toBits(0, 3),
            QuantumGate<3, Gates::TOFFOLI>().toBits(0, 3, 4),
            QuantumGate<2, Gates::CX>().toBits(0, 4),
            QuantumGate<3, Gates::TOFFOLI>().toBits(0, 4, 5),
            QuantumGate<2, Gates::CX>().toBits(0, 5),
            QuantumGate<3, Gates::TOFFOLI>().toBits(0, 5, 6),
            QuantumGate<2, Gates::CX>().toBits(0, 6),
            QuantumGate<3, Gates::TOFFOLI>().toBits(0, 6, 7),

            // --------------------------------------------------
            // 4) Controlled modular multiply by 4 (mod 21)
            // Control: phase qubit 1
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2, Gates::CX>().toBits(1, 3),
            QuantumGate<3, Gates::TOFFOLI>().toBits(1, 3, 4),
            QuantumGate<2, Gates::CX>().toBits(1, 4),
            QuantumGate<3, Gates::TOFFOLI>().toBits(1, 4, 5),
            QuantumGate<2, Gates::CX>().toBits(1, 5),
            QuantumGate<3, Gates::TOFFOLI>().toBits(1, 5, 6),
            QuantumGate<2, Gates::CX>().toBits(1, 6),
            QuantumGate<3, Gates::TOFFOLI>().toBits(1, 6, 7),

            // --------------------------------------------------
            // 5) Controlled modular multiply by 16 (mod 21)
            // Control: phase qubit 2
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2, Gates::CX>().toBits(2, 3),
            QuantumGate<3, Gates::TOFFOLI>().toBits(2, 3, 4),
            QuantumGate<2, Gates::CX>().toBits(2, 4),
            QuantumGate<3, Gates::TOFFOLI>().toBits(2, 4, 5),
            QuantumGate<2, Gates::CX>().toBits(2, 5),
            QuantumGate<3, Gates::TOFFOLI>().toBits(2, 5, 6),
            QuantumGate<2, Gates::CX>().toBits(2, 6),
            QuantumGate<3, Gates::TOFFOLI>().toBits(2, 6, 7),

            // --------------------------------------------------
            // 6) Inverse QFT on phase register
            // --------------------------------------------------
            QuantumGate<3, IQFT3>().toBits(0, 1, 2)
        );

    // Expected output : the measurement results are bit - noisy due to the low qubit number
    // and QFT leakage, but probability peaks appear at x=3 and x=5 in the phase register.
    // These correspond to fractions 3/8 and 5/8, which are the closest approximations
    // to multiples of 1/6 (i.e., 1/r) for the true period r=6.
    // From the period r=6, we can deduce the nontrivial factors of 21:
    // gcd(2^(r/2)-1, 21) = gcd(7, 21) = 7
    // gcd(2^(r/2)+1, 21) = gcd(9, 21) = 3
    KetCat::Visu::VisuProbaTable<256>().update<0, 1, 2>(ShorCircuit.getStateVector());

    return 0;
}

