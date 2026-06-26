#include "quantum_processor/qpu_diag.h"

int main()
{
    using namespace KetCat;
    using namespace SpectroscopicLetters;

    NeutralAtomTypeConfig
        <
        Element::Cs,

        256, /* Spatial discretization steps count */
        100.0, /* Spatial extent in a.u. */

        0, /* Index of the logical level 0 */
        2, /* Index of the logical level 1*/
        4, /* Index of the Rydberg level */

        QuantumNumber<6, s>,  /*0*/
        QuantumNumber<6, p>,  /*1*/
        QuantumNumber<7, s>,  /*2*/
        QuantumNumber<10, p>, /*3*/
        QuantumNumber<20, s>, /*4*/
        QuantumNumber<20, p>  /*5*/
        > Config;

    auto Grover = QuantumCircuit<2>().withGates(
		QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(1),

        QuantumGate<2, GateType::CZ>().toBits(0,1),

        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(1),

        QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::X>().toBits(1),
        QuantumGate<2, GateType::CZ>().toBits(0, 1),
        QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::X>().toBits(1),

        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(1)
    );
	
	auto Bell = QuantumCircuit<2>().withGates(
		QuantumGate<1, GateType::H>().toBits(0),
		QuantumGate<2, GateType::CX>().toBits(0, 1)
	);

    //QuantumProcessor<2, Config>("grover.kwf").execute(Grover); 
    QuantumProcessor<2, Config>("bell.kwf").execute(Bell);

    auto Diag = QPUDiagnostics<2, Config>::createQPUWithInitialState("test.kwf", 3);
    //Diag.QPU().execute(QuantumCircuit<2>().withGates(QuantumGate<2, GateType::CZ>().toBits(0, 1)));
}
