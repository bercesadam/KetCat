#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <array>

#include "quantum_processor/quantum_processor.h"
#include "gate_calibration/gate_diagnostics.h" 
#include "gate_calibration/phase_corrections.h"


namespace KetCat
{
    // Helper function to export the final static C++ lookup table matching the PolynomialInterpolator format
    template <natural_t SweepSteps>
    void generateCalibrationHeader(
        const std::array<real_t, SweepSteps>& sweepAngles,
        const std::array<real_t, SweepSteps>& rxErrors,
        const TwoQubitCalibResult& cz)
    {
        std::ofstream header("gate_calibration_map.h");
        header << "// ==========================================================================\n";
        header << "//  KetCat AUTOMATICALLY GENERATED OFFLINE CALIBRATION LOOKUP TABLE\n";
        header << "//  Execution date: 2026. \n";
        header << "// ==========================================================================\n";
        header << "#pragma once\n";
        header << "#include \"core_types.h\"\n";
        header << "#include \"gate_diagnostics/polynomial_interpolator.h\"\n\n";
        header << "namespace KetCat\n{\n";
        header << "    struct TwoQubitStaticCalib { real_t controlError; real_t targetError; real_t actualCzPhase; };\n\n";
        header << "    class GateCalibrationTable\n    {\n";
        header << "    public:\n";
        header << "        static constexpr natural_t CalibrationPoints = " << SweepSteps << ";\n\n";

        // Export Rx Interpolator Data
        header << "        // Exact Newton polynomial interpolator for arbitrary Rx rotations\n";
        header << "        static inline PolynomialInterpolator<CalibrationPoints> getRxCalib() noexcept\n";
        header << "        {\n";
        header << "            constexpr std::array<real_t, CalibrationPoints> x = {";
        for (size_t i = 0; i < SweepSteps; ++i) header << sweepAngles[i] << (i == SweepSteps - 1 ? "" : ", ");
        header << "};\n";
        header << "            constexpr std::array<real_t, CalibrationPoints> y = {";
        for (size_t i = 0; i < SweepSteps; ++i) header << rxErrors[i] << (i == SweepSteps - 1 ? "" : ", ");
        header << "};\n";
        header << "            return PolynomialInterpolator<CalibrationPoints>(x, y);\n";
        header << "        }\n\n";

        // Export CZ Data
        header << "        // Static corrections for the fixed two-qubit Rydberg CPhase (CZ) gate\n";
        header << "        static constexpr TwoQubitStaticCalib getCPhaseCalib() noexcept\n";
        header << "        {\n";
        header << "            return { " << cz.controlFramePhaseError << ", " << cz.targetFramePhaseError << ", " << cz.actualCzPhase << " };\n";
        header << "        }\n";
        header << "    };\n";
        header << "}\n";

        std::cout << "\n>> [KetCat] gate_calibration_map.h successfully generated!" << std::endl;
    }

    int main()
    {
        std::cout << "====================================================" << std::endl;
        std::cout << "     KetCat QUANTUM GATE SWEEP & NEWTON CALIBRATION " << std::endl;
        std::cout << "====================================================" << std::endl;

        using SpectroscopicLetters::s, SpectroscopicLetters::p;
        constexpr NeutralAtomTypeConfig
            <
            Element::Cs,

            256,  /* Spatial discretization steps count */
            100.0, /* Spatial extent in a.u. */

            0, /* Index of the logical level 0 */
            2, /* Index of the logical level 1 */
            5, /* Index of the Rydberg level */

            QuantumNumber<6, s>,  /*0*/
            QuantumNumber<6, p>,  /*1*/
            QuantumNumber<7, s>,  /*2*/
            QuantumNumber<10, p>, /*3*/
            QuantumNumber<20, s>, /*4*/
            QuantumNumber<20, p>  /*5*/
            > Config;

        constexpr natural_t sweep_steps = 10;
        std::array<real_t, sweep_steps> theta_points;
        std::array<real_t, sweep_steps> rx_phase_errors;
        std::array<real_t, sweep_steps> ry_phase_errors;

        // Prepare angles from 0 to pi using fixed-size containers required by the template
        for (size_t i = 0; i < sweep_steps; ++i)
        {
            real_t theta = (ConstexprMath::Pi * (i + 1)) / static_cast<real_t>(sweep_steps);
            theta_points[i] = theta;
        }

        using CircuitType = QuantumCircuit<1>;
        using QPUType = QuantumProcessor<1, Config>;
        using GlobalStateVectorType = decltype(std::declval<QPUType>().globalStateVector());

        // -------------------------------------------------------------------------
        // 1. SWEEP: Rx Gate calibration
        // -------------------------------------------------------------------------
        for (size_t i = 0; i < sweep_steps; ++i)
        {
            real_t theta = theta_points[i];
            std::array<GlobalStateVectorType, 2> basis_outputs;

            // Running from pure |0> state
            {
                auto Circuit = CircuitType().withGates(QuantumGate<1, GateType::RX>().withTheta(theta).toBits(0));
                auto QPU = QPUType("x_calib_0_" + std::to_string(theta) + ".kwf", 0);
                QPU.execute<1>(Circuit);
                basis_outputs[0] = QPU.globalStateVector();
            }

            // Running from pure |1> state
            {
                auto Circuit = CircuitType().withGates(QuantumGate<1, GateType::RX>().withTheta(theta).toBits(0));
                auto QPU = QPUType("x_calib_1_" + std::to_string(theta) + ".kwf", 1);
                QPU.execute<1>(Circuit);
                basis_outputs[1] = QPU.globalStateVector();
            }

            // Build the effective 2x2 matrix using the GateDiagnostic class
            Matrix<2> U_eff = GateDiagnostic<2>::buildEffectiveGate(basis_outputs);
            Matrix<2> U_ideal = GateDiagnostic<2>::buildIdealRx(theta);

            // Compute the phase error between the theoretical and effective representations
            real_t PhaseError = GateDiagnostic<2>::globalPhase(U_eff, U_ideal);
            rx_phase_errors[i] = PhaseError;

            std::cout << "Phase err at " << theta << " is: " << PhaseError << std::endl;
        }

        // -------------------------------------------------------------------------
        // HEADER EXPORT
        // -------------------------------------------------------------------------
        generateCalibrationHeader<sweep_steps>(theta_points, rx_phase_errors, TwoQubitCalibResult{});

        return 0;
    }
}

int main()
{
    KetCat::main();
    return 0;
}