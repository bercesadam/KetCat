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
    
    auto Circuit2 = QuantumCircuit<2>().withGates(
        QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<2, GateType::CZ>().toBits(0, 1),
        QuantumGate<2, GateType::CZ>().toBits(0, 1)
    );

    QuantumProcessor<2, Config>("smoke_test2.kwf").execute(Circuit2);

    auto Circuit3 = QuantumCircuit<2>().withGates(
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<2, GateType::CZ>().toBits(0, 1)
    );

    QuantumProcessor<2, Config>("smoke_test3.kwf").execute(Circuit3);


   
}