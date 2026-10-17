#include "config.h"
#include "quantum_processor/qpu_diag.h"

namespace KetCat::IonConfig
{
    using namespace SpectroscopicLetters;

    NeutralAtomTypeConfig
        <
        Element::Ca,

        456,   /* Spatial discretization steps count */
        750.0, /* Spatial extent in a.u. */

        0, /* Index of the logical level 0 */
        1, /* Index of the logical level 1*/
        2, /* currently unused */

        QuantumNumber<4, s>,  /* |g⟩ */
        QuantumNumber<3, d>,  /* |e⟩ */
        QuantumNumber<4, p>, 
        QuantumNumber<5, s>, 
        QuantumNumber<4, d>
        > Calcium_5Level;
}

int main()
{
    using namespace KetCat;

    auto Test = QuantumCircuit<1>().withGates(
        QuantumGate<1, GateType::Y>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::X>().toBits(0)
    );

	QPUDiagnostics<1, IonConfig::Calcium_5Level> DiagSession =
        QPUDiagnostics<1, IonConfig::Calcium_5Level>::createQPUWithInitialState("ion_test.kwf", std::bitset<1>{0});
	DiagSession.QPU().execute(Test);

    return(0);
}
