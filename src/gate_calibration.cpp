#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <array>

#define DISABLE_PHASE_GATE_CORRECTION
#include "config.h"
#include "quantum_processor/qpu_diag.h"
#include "gate_calibration/gate_diagnostics.h" 
#include "gate_calibration/phase_corrections.h"


namespace KetCat
{
    /// @brief Structure to hold static phase correction metadata for a two-qubit CPhase gate.
    struct TwoQubitStaticCalib
    {
        real_t controlError;    ///< Residual phase tracking error accumulated on the control qubit.
        real_t targetError;     ///< Residual phase tracking error accumulated on the target qubit.
        real_t actualCzPhase;   ///< The actual, physically realized entangling phase of the CZ gate.
    };

    /// @brief Generates an offline static C++ configuration header mapping physical calibration results.
    /// @details Outputs a structured header file (`gate_calibration_map.h`) containing a compiled-in 
    ///          Newton polynomial interpolator for arbitrary Rx single-qubit rotations and fixed corrections 
    ///          for the two-qubit Rydberg CPhase gate.
    /// @tparam SweepSteps The number of discrete sampling points evaluated during the calibration sweep.
    /// @param sweepAngles The collection of target rotation angles ($\theta$) used during the calibration sweep.
    /// @param rxErrors The calculated phase errors corresponding to each rotation angle in `sweepAngles`.
    /// @param cz The structured results of the two-qubit CZ entangling gate diagnostics.
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
        header << "#include \"gate_diagnostics/polynomial_interpolator.h\"\n\n#define CALIBRATED_GATES\n\n";
        header << "namespace KetCat\n{\n";
        header << "    class GateCalibrationTable\n    {\n";
        header << "    public:\n";
        header << "        static constexpr natural_t CalibrationPoints = " << SweepSteps + 1 << ";\n\n";

        // Export Rx Interpolator Data
        header << "        // Exact Newton polynomial interpolator for arbitrary Rx rotations\n";
        header << "        static inline PolynomialInterpolator<CalibrationPoints> getRxCalib() noexcept\n";
        header << "        {\n";
        header << "            constexpr std::array<real_t, CalibrationPoints> x = {0.0, ";
        for (size_t i = 0; i < SweepSteps; ++i) header << sweepAngles[i] << (i == SweepSteps - 1 ? "" : ", ");
        header << "};\n";
        header << "            constexpr std::array<real_t, CalibrationPoints> y = {0.0, ";
        for (size_t i = 0; i < SweepSteps; ++i) header << rxErrors[i] << (i == SweepSteps - 1 ? "" : ", ");
        header << "};\n";
        header << "            return PolynomialInterpolator<CalibrationPoints>(x, y);\n";
        header << "        }\n\n";

        // Export CZ Data
        header << "        // Static corrections for the fixed two-qubit Rydberg CPhase (CZ) gate\n";
        header << "        static constexpr TwoQubitStaticCalib getCPhaseCalib() noexcept\n";
        header << "        {\n";
        header << "            return { " << cz.m_controlFramePhaseError << ", " << cz.m_targetFramePhaseError << ", " << cz.m_actualCzPhase << " };\n";
        header << "        }\n";
        header << "    };\n";
        header << "}\n";

        std::cout << "\n>> [KetCat] gate_calibration_map.h successfully generated!" << std::endl;
    }

    /// @brief Executes a comprehensive calibration sweep for single-qubit Rx and two-qubit CZ gates.
    /// @details Simulates physical pulse actions on a neutral-atom hardware configuration using a multi-level 
    ///          atomic basis (e.g., Cesium). Quantifies phase deviations caused by off-resonant driving, 
    ///          Rydberg interactions, and Doppler shifts, exporting the results to an analytical lookup header.
    /// @return Returns 0 upon successful completion of the sweep sequence and header generation.
    int performGateCalibration()
    {
        auto Config = AtomConfig::Cesium_6Level;

        std::cout << "=========================================" << std::endl;
        std::cout << "     KetCat QUANTUM GATE CALIBRATION " << std::endl;
        std::cout << "=========================================" << std::endl;

        constexpr natural_t sweep_steps = 8;
        std::array<real_t, sweep_steps> ThetaSamplingPoints{};
        std::array<real_t, sweep_steps> RxGatePhaseErrors{};
        TwoQubitCalibResult CzResults{};

        // Prepare rotation angles ranging uniformly from 0 to Pi using fixed-size containers required by the template
        for (size_t i = 0; i < sweep_steps; ++i)
        {
            real_t theta = (ConstexprMath::Pi * (i + 1)) / static_cast<real_t>(sweep_steps);
            ThetaSamplingPoints[i] = theta;
        }

        using CircuitType = QuantumCircuit<1>;
        using QPUDiagType = QPUDiagnostics<1, Config>;
        using GlobalStateVectorType =
            StateVector<FiniteHilbertSpace<Config.LevelCount>, QuantumPicture::Schrodinger>;

        // -------------------------------------------------------------------------
        // 1. SWEEP: Single-Qubit Rx Gate Calibration
        // -------------------------------------------------------------------------
        for (size_t i = 0; i < sweep_steps; ++i)
        {
            real_t theta = ThetaSamplingPoints[i];
            std::array<GlobalStateVectorType, 2> BasisOutputs;
            auto Circuit = CircuitType().withGates(QuantumGate<1, GateType::RX>().withTheta(theta).toBits(0));

            // Execute the circuit using a pure |0> computational basis state initialization
            {
                auto Diag = QPUDiagType::createQPUWithInitialState("x_calib_0_" + std::to_string(theta) + ".kwf", 0);
                Diag.QPU().execute<1>(Circuit);
                BasisOutputs[0] = Diag.getGlobalStateVector();
                TimeMaster::Clock().reset();
            }

            // Execute the circuit using a pure |1> computational basis state initialization
            {
                auto Diag = QPUDiagType::createQPUWithInitialState("x_calib_1_" + std::to_string(theta) + ".kwf", 1);
                Diag.QPU().execute<1>(Circuit);
                BasisOutputs[1] = Diag.getGlobalStateVector();
                TimeMaster::Clock().reset();
            }

            // Reconstruct the empirical 2x2 unitary transformation matrix via the diagnostic subsystem
            Matrix<2> U_eff = GateDiagnostic<2>::buildEffectiveGate(BasisOutputs);
            Matrix<2> U_ideal = GateDiagnostic<2>::buildIdealRx(theta);

            // Compute the global/relative phase error discrepancy between the theoretical target and effective realizations
            real_t PhaseError = GateDiagnostic<2>::globalPhase(U_eff, U_ideal);
            RxGatePhaseErrors[i] = PhaseError;

            std::cout << "Phase err at " << theta << " is: " << PhaseError << std::endl;
        }

        // -------------------------------------------------------------------------
        // 2. SWEEP: Two-Qubit Controlled-Phase (CZ) Gate Diagnostics
        // -------------------------------------------------------------------------
        std::cout << "\n>> Starting Two-Qubit CZ Gate Phase-Tracking Sweep..." << std::endl;

        square_matrix_t<4> EffectiveCzGateMatrix{};

        // Iterate through all two-qubit computational basis states: |00>, |01>, |10>, and |11>
        for (natural_t BaseState = 0; BaseState < 4; ++BaseState)
        {
            auto Diag = QPUDiagnostics<2, Config>::createQPUWithInitialState("cz_calib_" + std::to_string(BaseState) + ".kwf", BaseState);
            Diag.QPU().execute(QuantumCircuit<2>().withGates(QuantumGate<2, GateType::CZ>().toBits(0, 1)));
            TimeMaster::Clock().reset();

            const auto EffectiveStateVector = Diag.getGlobalStateVector();

            // The total dimensionality of the 2-qubit global physical Hilbert space is 6 * 6 = 36 atomic levels.
            // We isolate and extract the complex amplitudes associated with the 4 pure logical computational basis states.
            // Using a standard little-endian bit sequencing layout, the physical index offsets are calculated below:
            constexpr natural_t Index00 = Config.Logical0Level + Config.Logical0Level * Config.LevelCount;
            constexpr natural_t Index01 = Config.Logical1Level + Config.Logical0Level * Config.LevelCount;
            constexpr natural_t Index10 = Config.Logical0Level + Config.Logical1Level * Config.LevelCount;
            constexpr natural_t Index11 = Config.Logical1Level + Config.Logical1Level * Config.LevelCount;

            EffectiveCzGateMatrix[0][BaseState] = EffectiveStateVector[Index00];
            EffectiveCzGateMatrix[1][BaseState] = EffectiveStateVector[Index01];
            EffectiveCzGateMatrix[2][BaseState] = EffectiveStateVector[Index10];
            EffectiveCzGateMatrix[3][BaseState] = EffectiveStateVector[Index11];
        }

        // Quantify the shifted azimuthal phases, cross-talk, and frame-tracking metrics
        CzResults = GateDiagnostic<4>::analyzeCPhase(EffectiveCzGateMatrix);

        // -------------------------------------------------------------------------
        // HEADER EXPORT
        // -------------------------------------------------------------------------
        generateCalibrationHeader<sweep_steps>(ThetaSamplingPoints, RxGatePhaseErrors, CzResults);

        return 0;
    }
}

int main()
{
    KetCat::performGateCalibration();
    return 0;
}
