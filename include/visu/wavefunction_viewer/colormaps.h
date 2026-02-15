#pragma once
#include <SDL2/SDL.h>
#include <algorithm>

namespace KetCat::Visu
{
    // Basic structure to represent an RGB color
    struct RGB
    {
        uint8_t r, g, b;
    };

    /**
     * Linearly interpolates between two RGB colors.
     * @param a The starting color.
     * @param b The ending color.
     * @param t The interpolation factor [0.0, 1.0].
     * @return The interpolated RGB color.
     */
    inline RGB lerp(RGB a, RGB b, double t)
    {
        return {
            static_cast<uint8_t>(a.r + (b.r - a.r) * t),
            static_cast<uint8_t>(a.g + (b.g - a.g) * t),
            static_cast<uint8_t>(a.b + (b.b - a.b) * t)
        };
    }

    /**
     * Inferno colormap: Perceptually uniform black-to-yellow heatmap.
     */
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

    /**
     * Viridis colormap: Popular perceptually uniform green-blue-yellow scale.
     */
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

    /**
     * Mako colormap: Dark blue to light teal, great for ocean/depth data.
     */
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

    /**
     * Flare colormap: Warm, glowing red-to-yellow transition.
     */
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

    /**
     * Vistia colormap: High-contrast yellow-orange-red to white scale.
     */
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
}
