#pragma once
#include <SDL2/SDL.h>
#include <algorithm>

///@file Color maps for SDl2 visualizations

namespace KetCat::Visu
{
    /// @brief Basic RGB color representation.
    struct RGB
    {
        uint8_t r, g, b;
    };

    /// @brief Linearly interpolate between two RGB colors.
    /// @param a Starting color.
    /// @param b Ending color.
    /// @param t Interpolation factor in [0.0, 1.0].
    /// @return The interpolated RGB color.
    inline RGB lerp(RGB a, RGB b, double t)
    {
        return {
            static_cast<uint8_t>(a.r + (b.r - a.r) * t),
            static_cast<uint8_t>(a.g + (b.g - a.g) * t),
            static_cast<uint8_t>(a.b + (b.b - a.b) * t)
        };
    }

    /// @brief Inferno colormap — perceptually uniform black→yellow heatmap.
    /// @param t Input value in [0, 1].
    /// @return Corresponding Inferno RGB color.
    inline RGB inferno(double t)
    {
        static const RGB c[] = {
            {0, 0, 4}, {31, 12, 72}, {85, 15, 109}, {136, 34, 106},
            {186, 54, 85}, {227, 89, 51}, {249, 140, 10}, {252, 255, 164}
        };

        t = std::clamp(t, 0.0, 1.0);
        double x = t * 7.0;
        int i = static_cast<int>(x);

        return lerp(c[i], c[std::min(i + 1, 7)], x - i);
    }

    /// @brief Viridis colormap — green→blue→yellow perceptually uniform scale.
    /// @param t Input value in [0, 1].
    /// @return Corresponding Viridis RGB color.
    inline RGB viridis(double t)
    {
        static const RGB c[] = {
            {68, 1, 84}, {59, 82, 139}, {33, 145, 140}, {94, 201, 98}, {253, 231, 37}
        };

        t = std::clamp(t, 0.0, 1.0);
        double x = t * 4.0;
        int i = static_cast<int>(x);

        return lerp(c[i], c[std::min(i + 1, 4)], x - i);
    }

    /// @brief Mako colormap — dark blue to light teal.
    /// @param t Input value in [0, 1].
    /// @return Corresponding Mako RGB color.
    inline RGB mako(double t)
    {
        static const RGB c[] = {
            {11, 4, 5}, {51, 46, 81}, {58, 111, 143}, {70, 185, 190}, {222, 252, 228}
        };

        t = std::clamp(t, 0.0, 1.0);
        double x = t * 4.0;
        int i = static_cast<int>(x);

        return lerp(c[i], c[std::min(i + 1, 4)], x - i);
    }

    /// @brief Flare colormap — warm red→orange→yellow gradient.
    /// @param t Input value in [0, 1].
    /// @return Corresponding Flare RGB color.
    inline RGB flare(double t)
    {
        static const RGB c[] = {
            {75, 25, 35}, {180, 75, 60}, {230, 145, 90}, {245, 210, 120}, {250, 245, 210}
        };

        t = std::clamp(t, 0.0, 1.0);
        double x = t * 4.0;
        int i = static_cast<int>(x);

        return lerp(c[i], c[std::min(i + 1, 4)], x - i);
    }

    /// @brief Vistia colormap — high‑contrast yellow→orange→red→white.
    /// @param t Input value in [0, 1].
    /// @return Corresponding Vistia RGB color.
    inline RGB vistia(double t)
    {
        static const RGB c[] = {
            {220, 255, 0}, {255, 190, 0}, {255, 125, 0}, {255, 60, 0}, {255, 255, 255}
        };

        t = std::clamp(t, 0.0, 1.0);
        double x = t * 4.0;
        int i = static_cast<int>(x);

        return lerp(c[i], c[std::min(i + 1, 4)], x - i);
    }

    /// @brief Scientific cyclic phase colormap.
    /// Red → Magenta → Blue → Cyan → Red (cyclic).
    /// @param phase Phase in [0,1] (cyclic).
    /// @param amplitude Amplitude in [0,1].
    inline RGB phase_colors(double phase)
    {
        phase = phase - std::floor(phase);       // cyclic wrap

        // Anchor colors
        const RGB red{ 255,   0,   0 };
        const RGB magenta{ 255,   0, 255 };
        const RGB blue{ 0,   0, 255 };
        const RGB cyan{ 0, 255, 255 };

        RGB color;

        if (phase < 0.25)
        {
            color = lerp(red, magenta, phase / 0.25);
        }
        else if (phase < 0.50)
        {
            color = lerp(magenta, blue, (phase - 0.25) / 0.25);
        }
        else if (phase < 0.75)
        {
            color = lerp(blue, cyan, (phase - 0.50) / 0.25);
        }
        else
        {
            color = lerp(cyan, red, (phase - 0.75) / 0.25);
        }

        return { color.r, color.g, color.b };
    }

    inline RGB phase_map(double phase, double amplitude)
    {
        amplitude = std::clamp(amplitude, 0.0, 1.0);
        RGB color = phase_colors(phase);

        // Amplitude as brightness
        return {
            static_cast<uint8_t>(color.r * amplitude),
            static_cast<uint8_t>(color.g * amplitude),
            static_cast<uint8_t>(color.b * amplitude)
        };
    }

}
