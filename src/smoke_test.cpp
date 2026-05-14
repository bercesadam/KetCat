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
        4, /* Index of the Rydberg level */

        QuantumNumber<6, s>,
        QuantumNumber<6, p>,
        QuantumNumber<7, s>,
        QuantumNumber<7, p>,
        QuantumNumber<40, s>,
        QuantumNumber<40, p>
        > Config;
        
    auto Circuit = QuantumCircuit<1>().withGates(
        QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::Y>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(0)
    );

    QuantumProcessor<1, Config>("smoke_test.kwf").execute(Circuit);
}