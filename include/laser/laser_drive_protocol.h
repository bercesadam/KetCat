#pragma once
#include <functional>
#include "atomic_units.h"
#include "hamiltonian/rabi_drive_hamiltonian.h"
#include "solvers/crank_nicolson_solver.h"
#include "systems/neutral_atom_manifold.h"


namespace KetCat
{
    /// @brief Time step for Time-Dependent Schrödinger Equation (TDSE) integration in atomic units (a.u.).
    static constexpr real_t TimeStepAu = 25;

    /// @brief Physical parameters for a laser pulse interaction.
    struct LaserPulse
    {
        real_t m_waveLengthNm;   ///< Wavelength in nanometers (nm).
        real_t m_intensityWCm2;  ///< Peak intensity in Watts per square centimeter (W/cm²).
    };

    /// @brief Quantization axes for Bloch sphere rotations.
    enum class RotationAxis
    {
        X, ///< Rotation around the X-axis (in-phase drive).
        Y, ///< Rotation around the Y-axis (quadrature drive, 90-degree shift).
        Z  ///< Rotation around the Z-axis (longitudinal phase shift/detuning).
    };

    /// @brief High-level command for a quantum gate operation.
    struct PulseCommand
    {
        RotationAxis m_axis;        ///< Target axis of rotation.
        real_t m_rotationAngleRad;  ///< Rotation angle θ in radians.
    };

    /// @brief Implements coherent control of a neutral atom using laser-driven pulses.
    ///
    /// @details
    /// This class handles the evolution of an atomic state vector under a multi-level
    /// Hamiltonian. It supports standard STIRAP-like pulse sequences for robust 
    /// population transfer and fractional pulses for arbitrary Bloch sphere rotations.
    ///
    /// @tparam Config Neutral atom configuration defining energy levels and dipole moments.
    template <NeutralAtomTypeConfig Config>
    class LaserDrivePulse
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using Manifold = NeutralAtomManifold<Config>;
        using HilbertSpace = typename Manifold::SingleAtomOperationHilbertSpace;

    public:
        /// @brief Callback for monitoring the state vector evolution during integration.
        using StepCallback = std::function<void(real_t time, const StateVector<HilbertSpace>&)>;

        /// @brief Simple Gaussian temporal profile for the Rabi frequency.
        ///
        /// @param t      Current time.
        /// @param t0     Peak position (mean).
        /// @param sigma  Width of the pulse (standard deviation).
        /// @return Pulse amplitude at time 't' relative to the peak.
        real_t gaussian(real_t t, real_t t0, real_t sigma)
        {
            real_t X = (t - t0) / sigma;
            return std::exp(-X * X);
        }

        /// @brief Executes a laser drive sequence on the target state vector.
        ///
        /// @details
        /// The function maps the `PulseCommand` to a specific STIRAP (Stimulated Raman 
        /// Adiabatic Passage) or detuning sequence. 
        /// - X/Y gates: Controlled via the relative phase between Pump and Stokes pulses.
        /// - Z gates: Implemented via a controlled energy shift (detuning) leading to 
        ///   dynamic phase accumulation.
        ///
        /// @param target      [in, out] The atomic state vector to evolve.
        /// @param baseLevel   Starting level index for the transition manifold.
        /// @param command     Rotation parameters (axis and angle).
        /// @param onStep      Optional callback for time-resolved state analysis.
        void operator()(StateVector<HilbertSpace>& target,
            const natural_t baseLevel,
            const PulseCommand& command,
            StepCallback onStep = nullptr)
        {
            const auto& Energies = Manifold::getHartreeEnergies();
            const auto& DipoleMatrix = Manifold::getDipoleMatrix();

            // Target Rabi frequency parameters for defining the pulse area.
            const real_t TargetRabiFrequency_Hz = 1e9;
            const real_t TargetRabiOmega_Au = Units::omegaAuFromHz(TargetRabiFrequency_Hz);
            const real_t Omega0 = TargetRabiOmega_Au;
            const real_t Sigma = 50.0 / Omega0;

            // Default simulation window for a full STIRAP sequence (approx 10σ).
            real_t FullTransitionTime_Au = 10.0 * Sigma;
            real_t PumpPhase = 0.0;
            real_t ZDetuning = 0.0;

            // --- Gate Parameter Logic ---
            if (command.m_axis == RotationAxis::X || command.m_axis == RotationAxis::Y)
            {
                // To rotate around Y, the Pump laser requires a π/2 phase shift 
                // relative to the Stokes reference.
                PumpPhase = (command.m_axis == RotationAxis::Y) ? (ConstexprMath::Pi / 2.0) : 0.0;

                // For arbitrary rotations (θ ≠ π), we use a "fractional STIRAP" approach.
                // We truncate the pulse sequence at the point where the mixing angle 
                // tan(θ/2) = Ω_pump / Ω_stokes matches the desired rotation.
                if (std::abs(command.m_rotationAngleRad - ConstexprMath::Pi) > 1e-6)
                {
                    real_t TargetTheta = command.m_rotationAngleRad / 2.0;

                    // Center-times for the delayed pulse sequence (Stokes precedes Pump).
                    real_t StokesTimeCenter = 4.0 * Sigma;
                    real_t PumpTimeCenter = 6.0 * Sigma;

                    // Calculate the truncation time T_limit by solving the ratio of two Gaussians.
                    // ln(tan(θ)) = ln(exp(-(t-tp)^2/s^2) / exp(-(t-ts)^2/s^2))
                    real_t DeltaTime = PumpTimeCenter - StokesTimeCenter;
                    FullTransitionTime_Au = StokesTimeCenter + (DeltaTime / 2.0) - (Sigma * Sigma / (2.0 * DeltaTime)) * std::log(std::tan(TargetTheta));
                }
            }
            else if (command.m_axis == RotationAxis::Z)
            {
                // Z-rotation: Pure phase accumulation.
                // Instead of driving transitions, we induce a temporary detuning Δ 
                // such that Δ * time = θ (rotation angle).
                ZDetuning = command.m_rotationAngleRad / (10.0 * Sigma);
            }

            const real_t t_stokes = 4.0 * Sigma;
            const real_t t_pump = 6.0 * Sigma;

            real_t Time = 0;
            natural_t TimeStep = 0;

            // --- TDSE Integration Loop ---
            while (Time < FullTransitionTime_Au)
            {
                // Calculate time-dependent Rabi frequencies.
                real_t Omega_stokes = Omega0 * gaussian(Time, t_stokes, Sigma);
                real_t Omega_pump = Omega0 * gaussian(Time, t_pump, Sigma);

                // For Z-rotations, disable pulse amplitudes to avoid population transfer.
                if (command.m_axis == RotationAxis::Z) { Omega_stokes = 0; Omega_pump = 0; }

                // Define the Rotating Wave Approximation (RWA) lasers.
                std::array<Laser, ConfigType::LevelCount - 1> lasers = { };

                // Laser 0 (Pump): Resonant with the first transition, phase-controlled.
                lasers[0] = { Energies[baseLevel + 1] - Energies[baseLevel],
                              Omega_pump / DipoleMatrix[0][1].re,
                              PumpPhase };

                // Laser 1 (Stokes): Resonant with the second transition + optional Z-detuning.
                lasers[1] = { (Energies[baseLevel + 2] - Energies[baseLevel + 1]) + ZDetuning,
                              Omega_stokes / DipoleMatrix[1][2].re,
                              0.0 };

                // Build the Hamiltonian and solve for the next time step using the 
                // unitary Crank-Nicolson method.
                MultiRwaRabiHamiltonian<ConfigType::LevelCount> H(Energies, DipoleMatrix, lasers);
                CrankNicolsonSolver<HilbertSpace> solver(H.getMatrix(), TimeStepAu);
                target = solver(target);

                Time += TimeStepAu;

                // Sparse callback invocation to reduce overhead.
                if (TimeStep % 1000000 == 0)
                {
                    if (onStep)
                    {
                        onStep(Time, target);
                    }
                }
                TimeStep++;
            }

            // Ensure the final state is captured.
            if (onStep) onStep(Time, target);
        }
    };
}
