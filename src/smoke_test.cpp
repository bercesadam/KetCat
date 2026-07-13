#include "config.h"
#include "quantum_processor/qpu_diag.h"


int main()
{
    using namespace KetCat;

    /// @brief 2-qubit Grover Search algorithm targeting the |11⟩ state.
    /// @details Applies uniform superposition, phase inversion via a Controlled-Z (CZ) gate, 
    /// followed by the standard diffusion operator (H, X, CZ, X, H).
    auto Grover = QuantumCircuit<2>().withGates(
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::H>().toBits(1),

        QuantumGate<2, GateType::CZ>().toBits(0, 1),

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
    

    /// @brief 3-qubit Greenberger-Horne-Zeilinger (GHZ) state preparation.
    /// @details Drives the neutral-atom registry into the maximally entangled state 
    /// |ψ⟩ = 1/√2 (|000⟩ + |111⟩) utilizing cascaded Controlled-NOT (CX) physical operations.
    auto GHZ = QuantumCircuit<3>().withGates(
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<2, GateType::CX>().toBits(0, 1),
        QuantumGate<2, GateType::CX>().toBits(1, 2)
    );


    /// @brief Evaluation of the rotation angle θ = 2 · acos(√(2/3)) for amplitude routing.
    constexpr real_t Theta_P0_2_3 = 2.0 * ConstexprMath::acos(ConstexprMath::sqrt(2.0 / 3.0));

    /// @brief 3-qubit Fair Dice state preparation circuit.
    /// @details Engineers a uniform superposition across a sub-manifold of six basis states, 
    /// mapping a random 6-sided dice roll. Employs asymmetric RY(θ) rotation angles and 
    /// entangled entangling networks to bound the outcome probabilities evenly.
    auto FairDice = QuantumCircuit<3>().withGates(
        QuantumGate<1, GateType::H>().toBits(0),
        QuantumGate<1, GateType::RX>().withTheta(Theta_P0_2_3).toBits(2),
        QuantumGate<1, GateType::X>().toBits(2),
        QuantumGate<1, GateType::RY>().withTheta(ConstexprMath::Pi / 4).toBits(1),
        QuantumGate<2, GateType::CX>().toBits(2, 1),
        QuantumGate<1, GateType::RY>().withTheta(7 * ConstexprMath::Pi / 4).toBits(1),
        QuantumGate<2, GateType::CX>().toBits(2, 1),
        QuantumGate<1, GateType::X>().toBits(2)
    );


    // Dispatch circuits onto the uncompressed Cesium 6-level physical hardware simulator engines
    QuantumProcessor<2, AtomConfig::Cesium_6Level>("grover.kwf").execute(Grover);
    QuantumProcessor<3, AtomConfig::Cesium_6Level>("ghz.kwf").execute(GHZ);
    QuantumProcessor<3, AtomConfig::Cesium_6Level>("dice.kwf").execute(FairDice);

    return(0);
}
