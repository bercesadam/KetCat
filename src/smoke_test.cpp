#include "config.h"
#include "quantum_processor/qpu_diag.h"


int main()
{
    using namespace KetCat;
    
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
	
	auto GHZ = QuantumCircuit<3>().withGates(
		QuantumGate<1, GateType::H>().toBits(0),
		QuantumGate<2, GateType::CX>().toBits(0, 1),
        QuantumGate<2, GateType::CX>().toBits(1, 2)
	);

    constexpr real_t Theta_P0_2_3 = 2.0 * ConstexprMath::acos(ConstexprMath::sqrt(2.0 / 3.0));
    auto FairDice = QuantumCircuit<3>().withGates(
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::RY>().withTheta(Theta_P0_2_3).toBits(2),
        QuantumGate<1, GateType::X>().toBits(2),
        QuantumGate<1, GateType::RY>().withTheta(ConstexprMath::Pi / 4).toBits(1),
        QuantumGate<2, GateType::CX>().toBits(2, 1),
        QuantumGate<1, GateType::RY>().withTheta(-ConstexprMath::Pi / 4).toBits(1),
        QuantumGate<2, GateType::CX>().toBits(2, 1),
        QuantumGate<1, GateType::X>().toBits(2)
    );

    auto Test = QuantumCircuit<1>().withGates(
        //QuantumGate<1, GateType::X>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(0)
    );

    //QuantumProcessor<1, AtomConfig::Cesium_6Level>("h.kwf").execute(Test);
    //QuantumProcessor<3, Config>("dice.kwf").execute(FairDice);
    //QuantumProcessor<2, AtomConfig::Cesium_6Level>("grover.kwf").execute(Grover);
    QuantumProcessor<3, AtomConfig::Cesium_6Level>("ghz.kwf").execute(GHZ);

    //auto Diag = QPUDiagnostics<2, Config>::createQPUWithInitialState("test.kwf", 3);
    //Diag.QPU().execute(QuantumCircuit<2>().withGates(QuantumGate<2, GateType::CZ>().toBits(0, 1)));
}
