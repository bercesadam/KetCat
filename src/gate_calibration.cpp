#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

#include "quantum_processor/quantum_processor.h"
#include "quantum_gates/gate_diag.h" 


namespace KetCat
{
    /// @brief Polynomial coefficients for quadratic fitting: y = c2*x^2 + c1*x + c0
    struct PolyCoefficients
    {
        real_t c2 = 0.0;
        real_t c1 = 0.0;
        real_t c0 = 0.0;
    };

    /// @brief Calibration result structure for a two-qubit gate
    struct TwoQubitCalibResult
    {
        real_t controlFramePhaseError = 0.0;
        real_t targetFramePhaseError = 0.0;
        real_t actualCzPhase = 0.0;
    };

    class GateCalibrator
    {
    public:
        /// @brief Analytical quadratic polynomial fitting (Least Squares Method)
        /// @details Numerically solves the X^T * X * C = X^T * Y system of normal equations using Cramer's rule.
        static PolyCoefficients fitQuadratic(const std::vector<real_t>& x, const std::vector<real_t>& y) noexcept
        {
            const size_t n = x.size();
            if (n < 3)
            {
                return PolyCoefficients{}; // A minimum of 3 points is required for a quadratic fit
            }

            real_t sum_x = 0.0, sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
            real_t sum_y = 0.0, sum_xy = 0.0, sum_x2y = 0.0;

            for (size_t i = 0; i < n; ++i)
            {
                real_t xi = x[i];
                real_t yi = y[i];
                real_t xi2 = xi * xi;

                sum_x += xi;
                sum_x2 += xi2;
                sum_x3 += xi2 * xi;
                sum_x4 += xi2 * xi2;
                sum_y += yi;
                sum_xy += xi * yi;
                sum_x2y += xi2 * yi;
            }

            // Determinant of the coefficient matrix (3x3)
            real_t det = n * (sum_x2 * sum_x4 - sum_x3 * sum_x3) -
                sum_x * (sum_x * sum_x4 - sum_x2 * sum_x3) +
                sum_x2 * (sum_x * sum_x3 - sum_x2 * sum_x2);

            if (std::abs(det) < 1E-12) return PolyCoefficients{}; // Numerically singular

            // Cramer's rule to extract the coefficients
            PolyCoefficients coeffs;

            // Determinant for c0 (constant term)
            real_t det_c0 = sum_y * (sum_x2 * sum_x4 - sum_x3 * sum_x3) -
                sum_x * (sum_xy * sum_x4 - sum_x2y * sum_x3) +
                sum_x2 * (sum_xy * sum_x3 - sum_x2y * sum_x2);

            // Determinant for c1 (linear term)
            real_t det_c1 = n * (sum_xy * sum_x4 - sum_x2y * sum_x3) -
                sum_y * (sum_x * sum_x4 - sum_x2 * sum_x3) +
                sum_x2 * (sum_x * sum_x2y - sum_x2 * sum_xy);

            // Determinant for c2 (quadratic term)
            real_t det_c2 = n * (sum_x2 * sum_x2y - sum_x3 * sum_xy) -
                sum_x * (sum_x * sum_x2y - sum_x2 * sum_xy) +
                sum_y * (sum_x * sum_x3 - sum_x2 * sum_x2);

            coeffs.c0 = det_c0 / det;
            coeffs.c1 = det_c1 / det;
            coeffs.c2 = det_c2 / det;

            return coeffs;
        }

        /// @brief Dynamically builds the theoretical, ideal unitary matrix of Rx(theta) matching GateTraits
        static Matrix<2> buildIdealRx(real_t theta) noexcept
        {
            const real_t halfTheta = theta / 2.0;
            Matrix<2> U;
            U.at(0, 0) = complex_t(std::cos(halfTheta), 0.0);
            U.at(0, 1) = complex_t(0.0, -std::sin(halfTheta));
            U.at(1, 0) = complex_t(0.0, -std::sin(halfTheta));
            U.at(1, 1) = complex_t(std::cos(halfTheta), 0.0);
            return U;
        }

        /// @brief Dynamically builds the theoretical, ideal unitary matrix of Ry(theta)
        static Matrix<2> buildIdealRy(real_t theta) noexcept
        {
            const real_t halfTheta = theta / 2.0;
            Matrix<2> U;
            U.at(0, 0) = complex_t(std::cos(halfTheta), 0.0);
            U.at(0, 1) = complex_t(-std::sin(halfTheta), 0.0);
            U.at(1, 0) = complex_t(std::sin(halfTheta), 0.0);
            U.at(1, 1) = complex_t(std::cos(halfTheta), 0.0);
            return U;
        }

        /// @brief Computes local phase errors and the pure CZ phase from the 4x4 effective matrix of CPhase
        static TwoQubitCalibResult analyzeCPhase(const square_matrix_t<4>& u_eff) noexcept
        {
            // Extract phases of the main diagonal: 0=|00>, 1=|01>, 2=|10>, 3=|11>
            real_t p00 = std::atan2(u_eff[0][0].im, u_eff[0][0].re);
            real_t p01 = std::atan2(u_eff[1][1].im, u_eff[1][1].re);
            real_t p10 = std::atan2(u_eff[2][2].im, u_eff[2][2].re);
            real_t p11 = std::atan2(u_eff[3][3].im, u_eff[3][3].re);

            TwoQubitCalibResult result;
            // Phase shift of the target qubit relative to |00>
            result.targetFramePhaseError = p01 - p00;
            // Phase shift of the control qubit relative to |00>
            result.controlFramePhaseError = p10 - p00;
            // The pure non-local two-qubit phase (ideally this is pi)
            result.actualCzPhase = p11 - p10 - p01 + p00;

            return result;
        }
    };


    // Helper function to export the final static C++ lookup table
    void generateCalibrationHeader(const PolyCoefficients& rx, const PolyCoefficients& ry, const TwoQubitCalibResult& cz)
    {
        std::ofstream header("gate_calibration_map.h");
        header << "// ==========================================================================\n";
        header << "//  KetCat AUTOMATICALLY GENERATED OFFLINE CALIBRATION LOOKUP TABLE\n";
        header << "//  Execution date: 2026. \n";
        header << "// ==========================================================================\n";
        header << "#pragma once\n";
        header << "#include \"core_types.h\"\n\n";
        header << "namespace KetCat\n{\n";
        header << "    struct PolyCoefficients { real_t c2; real_t c1; real_t c0; };\n";
        header << "    struct TwoQubitStaticCalib { real_t controlError; real_t targetError; real_t actualCzPhase; };\n\n";
        header << "    class GateCalibrationTable\n    {\n";
        header << "    public:\n";
        header << "        // Phase error polynomial coefficients for arbitrary Rx rotations\n";
        header << "        static constexpr PolyCoefficients getRxCalib() noexcept\n";
        header << "        {\n";
        header << "            return { " << rx.c2 << ", " << rx.c1 << ", " << rx.c0 << " };\n";
        header << "        }\n\n";
        header << "        // Phase error polynomial coefficients for arbitrary Ry rotations\n";
        header << "        static constexpr PolyCoefficients getRyCalib() noexcept\n";
        header << "        {\n";
        header << "            return { " << ry.c2 << ", " << ry.c1 << ", " << ry.c0 << " };\n";
        header << "        }\n\n";
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
        std::cout << "     KetCat QUANTUM GATE SWEEP & OLS CALIBRATION    " << std::endl;
        std::cout << "====================================================" << std::endl;

        using SpectroscopicLetters::s, SpectroscopicLetters::p;
        constexpr NeutralAtomTypeConfig
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

        constexpr size_t sweep_steps = 5;
        std::vector<real_t> theta_points;
        std::vector<real_t> rx_phase_errors;
        std::vector<real_t> ry_phase_errors;

        // Prepare angles from 0 to pi
        for (size_t i = 1; i <= sweep_steps; ++i)
        {
            real_t theta = (ConstexprMath::Pi * i) / static_cast<real_t>(sweep_steps);
            theta_points.push_back(theta);
        }

        using CircuitType = QuantumCircuit<1>;
        using QPUType = QuantumProcessor<1, Config>;
        using GlobalStateVectorType = decltype(std::declval<QPUType>().globalStateVector());

        // -------------------------------------------------------------------------
        // 1. SWEEP: Rx Gate calibration
        // -------------------------------------------------------------------------
        //std::cout << "\n[1/3] Running Rx gate sweep..." << std::endl;

        for (real_t theta : theta_points)
        {
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


            // Build the effective 2x2 matrix using your GateDiagnostic class
            Matrix<2> U_eff = GateDiagnostic<2>::buildEffectiveGate(basis_outputs);
            Matrix<2> U_ideal = GateCalibrator::buildIdealRx(theta);

            // Compute the phase error between the theoretical and effective representations
            real_t err = GateDiagnostic<2>::globalPhase(U_eff, U_ideal);
            rx_phase_errors.push_back(err);
        }

        for (real_t theta : theta_points)
        {
            std::array<GlobalStateVectorType, 2> basis_outputs;

            // Running from pure |1> state
            {
                auto Circuit = CircuitType().withGates(QuantumGate<1, GateType::RY>().withTheta(theta).toBits(0));
                auto QPU = QPUType("y_calib_1_" + std::to_string(theta) + ".kwf", 1);
                QPU.execute<1>(Circuit);
                basis_outputs[1] = QPU.globalStateVector();
            }

            // Running from pure |0> state
            {
                auto Circuit = CircuitType().withGates(QuantumGate<1, GateType::RY>().withTheta(theta).toBits(0));
                auto QPU = QPUType("y_calib_0_" + std::to_string(theta) + ".kwf", 0);
                QPU.execute<1>(Circuit);
                basis_outputs[0] = QPU.globalStateVector();
            }

            // Build the effective 2x2 matrix using your GateDiagnostic class
            Matrix<2> U_eff = GateDiagnostic<2>::buildEffectiveGate(basis_outputs);
            Matrix<2> U_ideal = GateCalibrator::buildIdealRy(theta);

            // Compute the phase error between the theoretical and effective representations
            real_t err = GateDiagnostic<2>::globalPhase(U_eff, U_ideal);
            ry_phase_errors.push_back(err);
        }

        PolyCoefficients rx_coeffs = GateCalibrator::fitQuadratic(theta_points, rx_phase_errors);
        std::cout << " -> Rx Polynomial: c2=" << rx_coeffs.c2 << ", c1=" << rx_coeffs.c1 << ", c0=" << rx_coeffs.c0 << std::endl;

        PolyCoefficients ry_coeffs = GateCalibrator::fitQuadratic(theta_points, ry_phase_errors);
        std::cout << " -> Ry Polynomial: c2=" << ry_coeffs.c2 << ", c1=" << ry_coeffs.c1 << ", c0=" << ry_coeffs.c0 << std::endl;
        /*
        // -------------------------------------------------------------------------
        // 3. CALIBRATION: Two-qubit CPhase (36x36 Hamiltonian space)
        // -------------------------------------------------------------------------
        std::cout << "\n[3/3] Calibrating CPhase (36x36 Hamiltonian) two-atom gate..." << std::endl;
        std::array<StateVector<TwoAtom36LevelSpace>, 4> basis_outputs_2q;
        std::array<natural_t, 4> logical_indices = { 0, 2, 12, 14 }; // Computational logical basis indices in your 36-level space

        for (size_t col = 0; col < 4; ++col)
        {
            QpuSimulator36x36 sim2q;
            sim2q.setInitialState(logical_indices[col]);
            sim2q.runCPhasePulseSequence(); // STIRAP -> Blockade -> Inverted STIRAP
            basis_outputs_2q[col] = sim2q.getGlobalStateVector();
        }

        // Extract the 4x4 effective logical matrix using your diagnostic utilities
        square_matrix_t<4> U_eff_2q = GateDiagnostic<4>::buildEffectiveGate(basis_outputs_2q);
        TwoQubitCalibResult cz_result = GateCalibrator::analyzeCPhase(U_eff_2q);

        std::cout << " -> Control local error: " << cz_result.controlFramePhaseError << " rad." << std::endl;
        std::cout << " -> Target local error: " << cz_result.targetFramePhaseError << " rad." << std::endl;
        std::cout << " -> Actual Rydberg phase twist: " << cz_result.actualCzPhase << " rad." << std::endl;
        */
        // -------------------------------------------------------------------------
        // HEADER EXPORT
        // -------------------------------------------------------------------------
        generateCalibrationHeader(rx_coeffs, ry_coeffs, TwoQubitCalibResult{});
        
        return 0;
    }
}

int main()
{
    KetCat::main();
    return 0;
}