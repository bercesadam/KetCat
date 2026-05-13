#include "quantum_processor/neutral_atom_computer.h"

int main()
{
    using namespace KetCat;
    using namespace SpectroscopicLetters;

    real_t TimeStep = 15;

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
        
    auto qc = QuantumCircuit<1>().withGates(
        QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::Y>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::RY>().withTheta(ConstexprMath::Pi / 4).toBits(0)
    );
}