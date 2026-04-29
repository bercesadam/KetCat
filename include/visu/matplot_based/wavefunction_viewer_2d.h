#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <matplot/matplot.h>

#include "core_types.h"
#include "hilbert_space/hilbert.h"
#include "hilbert_space/state_vector.h"

namespace KetCat::Visu
{
    /// @brief Matplot++ alapú vizualizáló 2D hullámfüggvényekhez
    template<typename HilbertSpace>
    class WavefunctionViewer2D
    {
    private:
        static constexpr int Grid = HilbertSpace::Steps;

        matplot::figure_handle m_Fig;
        matplot::axes_handle m_AxProb, m_AxReal, m_AxImag, m_AxPhase;

        bool m_Initialized = false;


        std::vector<std::vector<double>> x_coords;
        std::vector<std::vector<double>> y_coords;
        std::vector<std::vector<double>> prob;
        std::vector<std::vector<double>> re;
        std::vector<std::vector<double>> im;
        std::vector<std::vector<double>> phase;


        matplot::contours_handle  probPlot;
        matplot::contours_handle  rePlot;
        matplot::contours_handle  imPlot;
        matplot::contours_handle  phasePlot;

    public:
        WavefunctionViewer2D(int width = 1200, int height = 800)
        {
            // Környezeti változók és inicializálás (mint az 1D-s kódban)
#ifdef _WIN32
            _putenv_s("GNUTERM", "png");
#else
            setenv("GNUTERM", "png", 1);
#endif

            m_Fig = matplot::figure(false);
            m_Fig->size(width, height);

            // 2x2-es elrendezés létrehozása
            m_AxProb = m_Fig->add_subplot(2, 2, 0);
            m_AxImag = m_Fig->add_subplot(2, 2, 1);
            m_AxPhase = m_Fig->add_subplot(2, 2, 2);
            m_AxReal = m_Fig->add_subplot(2, 2, 3);

            // 1. Koordináta tengelyek létrehozása (X: 0..Grid-1, Y: 0..Grid-1)
            std::vector<double> x_coords = matplot::linspace(0, Grid - 1, Grid);
            std::vector<double> y_coords = matplot::linspace(0, Grid - 1, Grid);

            // 2. 1D-be kiterített (flat) vektorok a Z adatoknak
            prob = std::vector<std::vector<double>>(Grid, std::vector<double>(Grid));
            re = std::vector<std::vector<double>>(Grid, std::vector<double>(Grid));
            im = std::vector<std::vector<double>>(Grid, std::vector<double>(Grid));
            phase = std::vector<std::vector<double>>(Grid, std::vector<double>(Grid));

        }

        template<typename StateVectorType>
        void render(const StateVectorType& psi, const std::string& customTitle)
        {
            for (int y = 0; y < Grid; ++y) {
                for (int x = 0; x < Grid; ++x) {
                    auto val = psi[{static_cast<natural_t>(x), static_cast<natural_t>(y)}];

                    // Matplot++ mátrix indexelés: [sor][oszlop] -> [y][x]
                    prob[y][x] = val.normSquared();
                    re[y][x] = val.re;
                    im[y][x] = val.im;
                    phase[y][x] = std::atan2(val.im, val.re);
                }
            }

            m_Fig->title(customTitle);

            // Első futásnál beállítjuk a hőtérképeket és a színskálákat
            if (!m_Initialized) {
                probPlot = setupAxes(m_AxProb, prob, "Probability |ψ|²", matplot::palette::inferno());
                imPlot = setupAxes(m_AxImag, im, "Im(ψ)", matplot::palette::summer());
                phasePlot = setupAxes(m_AxPhase, phase, "Phase arg(ψ)", matplot::palette::hsv());
                rePlot = setupAxes(m_AxReal, re, "Re(ψ)", matplot::palette::winter());

                m_Initialized = true;
            }
            else {
                updateAxes(probPlot, prob);
                updateAxes(imPlot, im);
                updateAxes(phasePlot, phase);
                updateAxes(rePlot, re);
            }

            static int frameCounter = 0;
            std::string filename = "frame_2d_" + std::to_string(frameCounter++) + ".png";
            m_Fig->save(filename);
        }

    private:
        matplot::contours_handle setupAxes(matplot::axes_handle ax, const std::vector<std::vector<double>>& data,
            const std::string& label, const matplot::vector_2d& cmap)
        {
            ax->hold(matplot::off);
            auto hnd = ax->contourf(x_coords, y_coords, data);
            ax->title(label);
            ax->colormap(cmap);
            ax->axis(matplot::square);
            return hnd;
        }

        void updateAxes(matplot::contours_handle hnd, const std::vector <std::vector<double>>& data)
        {
            hnd->z_data(data);
        }
    };
}
