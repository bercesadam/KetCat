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

    auto Circuit = QuantumCircuit<2>().withGates(
		QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<2, GateType::CX>().toBits(0,1)
    );

	//QuantumProcessor<2, Config>("smoke_test_1").execute(Circuit);

    auto DiagSession = QPUDiagnostics<2, Config>::createQPUWithInitialState("smoke_test", std::bitset<2>{"00"});
    DiagSession.QPU().execute(Circuit);
  
}