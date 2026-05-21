#include "quantum_processor/quantum_processor.h"

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
        5, /* Index of the Rydberg level */

        QuantumNumber<6, s>,  /*0*/
        QuantumNumber<6, p>,  /*1*/
        QuantumNumber<7, s>,  /*2*/
        QuantumNumber<10, p>, /*3*/
        QuantumNumber<20, s>, /*4*/
        QuantumNumber<20, p>  /*5*/
        > Config;
        
    auto Circuit = QuantumCircuit<1>().withGates(
        QuantumGate<1, GateType::X>().toBits(0),
		QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::Y>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(0)
    );

    QuantumProcessor<1, Config>("smoke_test.kwf").execute(Circuit);
}