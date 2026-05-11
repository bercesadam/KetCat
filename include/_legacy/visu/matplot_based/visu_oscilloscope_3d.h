#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <string>
#include <iostream>

#include "matplot/matplot.h"
#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "hilbert_space/state_vector.h"

namespace KetCat::Visu
{

    /// @brief Helper function: set an environment variable in a cross-platform way
    inline void set_env(const std::string& name, const std::string& value) {
    #ifdef _WIN32
        _putenv_s(name.c_str(), value.c_str());
    #else
        setenv(name.c_str(), value.c_str(), 1);
    #endif
    }

    /// @brief High-fidelity 3D visualization for 1D quantum states using Matplotplusplus
    ///
    /// @details
    /// The visualization includes:
    ///  - A 3D complex plane view (Re(ψ) and Im(ψ) vs Position)
    ///  - Probability density plot |ψ|²
    ///  - Phase angle plot in radians
    ///  - Status overlay for simulation metadata
    class VisuOscilloscope3D
    {
    private:
        matplot::figure_handle m_Fig;
        matplot::labels_handle m_StatusLabel;

        matplot::axes_handle m_Ax3d;
        matplot::line_handle m_WavefunctionLine;
        matplot::line_handle m_AxisLine;

        matplot::axes_handle m_AxProb;
        matplot::line_handle m_ProbLine;

        matplot::axes_handle m_AxPhase;
        matplot::line_handle m_PhaseLine;

        bool m_Initialized = false;

    public:
        /// @brief Construct a VisuOscilloscope3D and initialize the plotting layout
        VisuOscilloscope3D()
        {
            set_env("GNUTERM", "png");
            m_Fig = matplot::figure(false);
            m_Fig->backend()->run_command("set terminal pngcairo size 1200,800");
            m_Fig->backend()->run_command("set output");
            m_Fig->size(1200, 800);

            m_Ax3d = m_Fig->add_subplot(3, 2, { 0, 1, 2, 3 });
            m_Ax3d->view(45, 30);
            m_Ax3d->xlabel("Position (x)");
            m_Ax3d->ylabel("Re(ψ)");
            m_Ax3d->zlabel("Im(ψ)");
            m_Ax3d->grid(true);

            m_AxProb = m_Fig->add_subplot(3, 2, 4);
            m_AxProb->title("Probability Density |ψ|²");
            m_AxProb->grid(true);

            m_AxPhase = m_Fig->add_subplot(3, 2, 5);
            m_AxPhase->title("Phase Angle (rad)");

            double Pi = 3.14159265358979;
            m_AxPhase->yticks({ -Pi, -Pi / 2, 0, Pi / 2, Pi });
            m_AxPhase->yticklabels({ "-π", "-π/2", "0", "π/2", "π" });
            m_AxPhase->ylim({ -3.2, 3.2 });
            m_AxPhase->grid(true);

            auto AxOverlay = m_Fig->add_axes({ 0, 0, 1, 1 });
            AxOverlay->xlim({ 0, 1 });
            AxOverlay->ylim({ 0, 1 });
            AxOverlay->axis(false);
            AxOverlay->hold(matplot::on);

            m_StatusLabel = AxOverlay->text(0.02, 0.96, "Initializing...");
            m_StatusLabel->font_size(14);
            m_StatusLabel->color("black");

            m_Initialized = false;
        }

        /// @brief Update the visualization with a new state vector and save as frame
        /// @tparam HilbertSpace Type of the spatial Hilbert space (must be 1D)
        /// @param psi The current state vector to visualize
        /// @param statusText Metadata string to display on the overlay
        template<KetCat::spatial_hilbert_space_with_dim_t<1_D> HilbertSpace>
        void update(const StateVector<HilbertSpace>& psi, const std::string& statusText) {
            size_t N = HilbertSpace::Dim;

            std::vector<double> X(N), Re(N), Im(N), Prob(N), Phase(N);
            for (size_t I = 0; I < N; ++I) {
                X[I] = static_cast<double>(I) / (N - 1);
                Re[I] = psi[I].re;
                Im[I] = psi[I].im;
                Prob[I] = psi[I].normSquared();
                Phase[I] = std::atan2(Im[I], Re[I]);
            }

            m_StatusLabel->label_values({ statusText });

            if (!m_Initialized) {
                m_Ax3d->hold(matplot::on);
                m_WavefunctionLine = m_Ax3d->plot3(X, Re, Im, "b-");
                m_WavefunctionLine->line_width(2);

                std::vector<double> AxisX = { 0, 1 };
                std::vector<double> AxisZeros = { 0, 0 };
                m_AxisLine = m_Ax3d->plot3(AxisX, AxisZeros, AxisZeros, "r-");
                m_AxisLine->line_width(4);
                m_Ax3d->hold(matplot::off);

                m_ProbLine = m_AxProb->plot(X, Prob, "c-");
                m_ProbLine->line_width(2);

                m_PhaseLine = m_AxPhase->plot(X, Phase, "g-");
                m_PhaseLine->line_width(2);

                m_Initialized = true;
            }
            else {
                m_WavefunctionLine->x_data(X);
                m_WavefunctionLine->y_data(Re);
                m_WavefunctionLine->z_data(Im);
                m_ProbLine->y_data(Prob);
                m_PhaseLine->y_data(Phase);
            }

            static int FrameCounter = 0;
            std::string Filename = "frame_" + std::to_string(FrameCounter++) + ".png";
            m_Fig->save(Filename);
        }

        /// @brief Set the camera view angle for the 3D complex plane
        /// @param azimuth Azimuthal angle in degrees
        /// @param elevation Elevation angle in degrees
        void setView(double azimuth, double elevation) {
            m_Ax3d->view(azimuth, elevation);
        }

        /// @brief Set the axis limits for the real and imaginary parts
        /// @param reMin Minimum value for the real axis
        /// @param reMax Maximum value for the real axis
        /// @param imMin Minimum value for the imaginary axis
        /// @param imMax Maximum value for the imaginary axis
        void setLimits(double reMin, double reMax, double imMin, double imMax) {
            m_Ax3d->ylim({ reMin, reMax });
            m_Ax3d->zlim({ imMin, imMax });
        }
    };
} // namespace KetCat::Visu
