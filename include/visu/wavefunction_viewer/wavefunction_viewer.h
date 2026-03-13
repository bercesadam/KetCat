#pragma once

#include <SDL.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_image.h>

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <iomanip>

#include "colormaps.h"
#include "hilbert_space/hilbert.h"
#include "hilbert_space/state_vector.h"

namespace KetCat::Visu
{
    /// @brief SDL‑based visualizer for 2D hydrogen wavefunctions.
    /// @details
    /// Displays four panels:
    ///  • Probability density |ψ|²  
    ///  • Imaginary part Im(ψ)  
    ///  • Phase arg(ψ)  
    ///  • Real part Re(ψ)  
    ///
    /// Each panel uses its own colormap and color bar.
    ///
    /// @tparam HilbertSpace  
    /// Must provide:  
    ///  • static constexpr  int Steps  
    ///  • static  int getIndex( int x,  int y)
    template<typename HilbertSpace>
    class WavefunctionViewer
    {
        static constexpr int Grid = HilbertSpace::Steps;

        /// @brief Window width in pixels
         int m_width;

        /// @brief Window height in pixels
         int m_height;

        /// @brief Frame counter used for automatic PNG saving
         int m_frameCounter;

        /// @brief SDL window and renderer handles
        SDL_Window* m_window{};
        SDL_Renderer* m_renderer{};

        /// @brief Fonts for text rendering
        TTF_Font* m_font{};
        TTF_Font* m_smallFont{};

    public:
        /// @brief Construct the visualization window.
        /// @param width  Window width in pixels.
        /// @param height Window height in pixels.
        WavefunctionViewer( int width,  int height)
            : m_width(width), m_height(height), m_frameCounter(0)
        {
            SDL_Init(SDL_INIT_VIDEO);
            TTF_Init();

            m_font = TTF_OpenFont(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18);

            m_smallFont = TTF_OpenFont(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14);

            m_window = SDL_CreateWindow(
                "Wavefunction Viewer",
                SDL_WINDOWPOS_CENTERED,
                SDL_WINDOWPOS_CENTERED,
                width, height, 0);

            m_renderer = SDL_CreateRenderer(
                m_window, -1, SDL_RENDERER_ACCELERATED);
        }

        /// @brief Destructor: releases SDL, TTF and window resources.
        ~WavefunctionViewer()
        {
            TTF_CloseFont(m_font);
            TTF_CloseFont(m_smallFont);
            TTF_Quit();

            SDL_DestroyRenderer(m_renderer);
            SDL_DestroyWindow(m_window);
            SDL_Quit();
        }

        /// @brief Render all wavefunction panels and overlay information.
        ///
        /// @tparam StateVectorType  
        /// Must provide:  
        ///  • operator[]( int)  
        ///  • members .re, .im, .normSquared()
        ///
        /// @param psi   The wavefunction ψ.
        /// @param alpha Complex coefficient of state |ψ₀⟩.
        /// @param beta  Complex coefficient of state |ψ₁⟩.
        template<typename StateVectorType>
        void render(const StateVectorType& psi,
                    std::string customTitle)
        {
            SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
            SDL_RenderClear(m_renderer);

            const  int margin = 40;
            const  int barSpace = 50;
            const  int boxW =
                (m_width - (margin * 3) - (barSpace * 2)) / 2;
            const  int boxH =
                (m_height - margin * 3) / 2;

            // --- Density |ψ|² (top-left)
             int dX = margin;
             int dY = margin;
            drawDensity(psi, dX, dY, boxW, boxH);
            drawColorBar(
                dX + boxW + 5, dY, 15, boxH,
                inferno,
                { {1.0, "max"}, {0.5, "0.5"}, {0.0, "0"} });

            // --- Imag part Im(ψ) (top-right)
             int iX = margin * 2 + boxW + barSpace;
             int iY = margin;
            drawImag(psi, iX, iY, boxW, boxH);
            drawColorBar(
                iX + boxW + 5, iY, 15, boxH,
                flare,
                { {1.0, "+max"}, {0.5, "0"}, {0.0, "-max"} });

            // --- Phase arg(ψ) (bottom-left)
             int pX = margin;
             int pY = margin * 2 + boxH;
            drawPhase(psi, pX, pY, boxW, boxH);
            drawColorBar(
                pX + boxW + 5, pY, 15, boxH,
                phase_colors,
                { {1.0,  "π"},
                  {0.75, "π/2"},
                  {0.5,  "0"},
                  {0.25, "-π/2"},
                  {0.0,  "-π"} });

            // --- Real part Re(ψ) (bottom-right)
             int rX = margin * 2 + boxW + barSpace;
             int rY = margin * 2 + boxH;
            drawReal(psi, rX, rY, boxW, boxH);
            drawColorBar(
                rX + boxW + 5, rY, 15, boxH,
                mako,
                { {1.0, "+max"}, {0.5, "0"}, {0.0, "-max"} });

            // Render panel labels
            drawLabels(margin, boxW, boxH, barSpace);

            // Render externally set custom title
            drawText(customTitle, 1, 1,
                     SDL_Color{255,255,255,255},
                     m_font);

            SDL_RenderPresent(m_renderer);

            ++m_frameCounter;
            saveFrame();
        }

    private:
        /// @brief Draw UTF‑8 text using SDL_ttf.
        void drawText(const std::string& text,
                       int x,
                       int y,
                      SDL_Color color,
                      TTF_Font* font)
        {
            SDL_Surface* surf =
                TTF_RenderUTF8_Blended(font, text.c_str(), color);

            SDL_Texture* tex =
                SDL_CreateTextureFromSurface(m_renderer, surf);

            SDL_Rect dst{ x, y, surf->w, surf->h };
            SDL_RenderCopy(m_renderer, tex, nullptr, &dst);

            SDL_FreeSurface(surf);
            SDL_DestroyTexture(tex);
        }

        /// @brief Draw panel labels below each visualization box.
        void drawLabels( int margin,
                         int boxW,
                         int boxH,
                         int barSpace)
        {
            SDL_Color white{255,255,255,255};
             int col2X = margin * 2 + boxW + barSpace;

            drawText("|ψ|²",
                     margin,
                     margin + boxH + 5,
                     white,
                     m_font);

            drawText("Im(ψ)",
                     col2X,
                     margin + boxH + 5,
                     white,
                     m_font);

            drawText("arg(ψ)",
                     margin,
                     margin * 2 + boxH + boxH + 5,
                     white,
                     m_font);

            drawText("Re(ψ)",
                     col2X,
                     margin * 2 + boxH + boxH + 5,
                     white,
                     m_font);
        }

        /// @brief Draw one pixel.
        void setPixel( int x,  int y, RGB c)
        {
            SDL_SetRenderDrawColor(m_renderer,
                                   c.r, c.g, c.b, 255);
            SDL_RenderDrawPoint(m_renderer, x, y);
        }

        /// @brief Draw probability density |ψ|².
        template<typename StateVectorType>
        void drawDensity(const StateVectorType& psi,
                          int ox,  int oy,
                          int w,  int h)
        {
            double maxVal = 0.0;

            for (int i = 0; i < Grid * Grid; ++i)
                maxVal = std::max(maxVal, psi[i].normSquared());

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    natural_t gx = x * Grid / w;
                    natural_t gy = y * Grid / h;

                    double v = psi[{gx, gy}].normSquared() / maxVal;

                    RGB c = inferno(std::tanh(3.0 * v));
                    setPixel(ox + x, oy + y, c);
                }
            }
        }

        /// @brief Draw real part Re(ψ).
        template<typename StateVectorType>
        void drawReal(const StateVectorType& psi,
                       int ox,  int oy,
                       int w,  int h)
        {
            double maxAbs = 0.0;

            for (int i = 0; i < Grid * Grid; ++i)
                maxAbs = std::max(maxAbs, std::abs(psi[i].re));

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    natural_t gx = x * Grid / w;
                    natural_t gy = y * Grid / h;
                    double v = psi[{gx, gy}].re / (maxAbs + 1e-12);

                    double val = 0.5 + 0.5 * std::tanh(2.5 * v);

                    RGB c = mako(val);
                    setPixel(ox + x, oy + y, c);
                }
            }
        }

        /// @brief Draw imaginary part Im(ψ).
        template<typename StateVectorType>
        void drawImag(const StateVectorType& psi,
                       int ox,  int oy,
                       int w,  int h)
        {
            double maxAbs = 0.0;

            for (int i = 0; i < Grid * Grid; ++i)
                maxAbs = std::max(maxAbs, std::abs(psi[i].im));

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    natural_t gx = x * Grid / w;
                    natural_t gy = y * Grid / h;

                    double v = psi[{gx, gy}].im / (maxAbs + 1e-12);

                    double val = 0.5 + 0.5 * std::tanh(2.5 * v);

                    RGB c = flare(val);
                    setPixel(ox + x, oy + y, c);
                }
            }
        }

        /// @brief Draw phase arg(ψ) mapped from [−π, π] to [0, 1].
        template<typename StateVectorType>
        void drawPhase(const StateVectorType& psi,
                        int ox,  int oy,
                        int w,  int h)
        {
            double maxVal = 0.0;

            for (int i = 0; i < Grid * Grid; ++i)
                maxVal = std::max(maxVal, psi[i].normSquared());

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    natural_t gx = x * Grid / w;
                    natural_t gy = y * Grid / h;

                    complex_t a = psi[{gx, gy}];

                    double ampl = a.normSquared() / maxVal;

                    double angle =
                        std::atan2(a.im, a.re);

                    double phase = (angle + M_PI) / (2.0 * M_PI);

                    RGB c = phase_map(phase, ampl);
                    setPixel(ox + x, oy + y, c);
                }
            }
        }

        /// @brief Draw a vertical color bar with labels.
        /// @param ox        X position.
        /// @param oy        Y position.
        /// @param w         Width in pixels.
        /// @param h         Height in pixels.
        /// @param colormap  Function mapping [0,1] → RGB.
        /// @param labels    Vector of (relative position, text).
        void drawColorBar(
             int ox,  int oy,
             int w,  int h,
            RGB(*colormap)(double),
            const std::vector<std::pair<double, std::string>>& labels)
        {
            for ( int y = 0; y < h; ++y)
            {
                double t = 1.0 - (static_cast<double>(y) / h);
                RGB c = colormap(t);

                SDL_SetRenderDrawColor(
                    m_renderer,
                    c.r, c.g, c.b, 255);

                SDL_RenderDrawLine(
                    m_renderer,
                    ox,     oy + y,
                    ox + w, oy + y);
            }

            SDL_SetRenderDrawColor(
                m_renderer, 100,100,100,255);

            SDL_Rect border{ ox, oy, w, h };
            SDL_RenderDrawRect(m_renderer, &border);

            SDL_Color white{255,255,255,255};

            for (const auto& label : labels)
            {
                 int labelY =
                    oy + static_cast< int>(
                        (1.0 - label.first) * h) - 7;

                drawText(label.second,
                         ox + w + 8,
                         labelY,
                         white,
                         m_smallFont);

                SDL_RenderDrawLine(
                    m_renderer,
                    ox + w,
                    labelY + 7,
                    ox + w + 4,
                    labelY + 7);
            }
        }

        /// @brief Save the current frame as a PNG image.
        void saveFrame()
        {
            SDL_Surface* surface =
                SDL_CreateRGBSurfaceWithFormat(
                    0,
                    m_width,
                    m_height,
                    32,
                    SDL_PIXELFORMAT_RGBA32);

            if (!surface)
                return;

            SDL_RenderReadPixels(
                m_renderer,
                nullptr,
                SDL_PIXELFORMAT_RGBA32,
                surface->pixels,
                surface->pitch);

            std::ostringstream filename;
            filename << std::setw(4)
                     << std::setfill('0')
                     << m_frameCounter
                     << "_.png";

            IMG_SavePNG(surface, filename.str().c_str());

            SDL_FreeSurface(surface);
        }
    };

} // namespace KetCat::Visu
