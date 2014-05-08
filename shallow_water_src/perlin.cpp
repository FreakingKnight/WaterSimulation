/*
 * FILE:
 *   perlin.cpp
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * CREATED:
 *   25-Oct-2011
 *
 * COPYRIGHT:
 *   Copyright (C) 2012, Stephen Thompson.
 *
 *   This file is part of Stephen Thompson's Shallow Water Demo.
 *
 *   The Shallow Water Demo is free software: you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation, either
 *   version 3 of the License, or (at your option) any later version.
 *
 *   The Shallow Water Demo is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the Shallow Water Demo. If not, see
 *   <http://www.gnu.org/licenses/>.
 *
 */

#include "perlin.hpp"

#include <cmath>
#include <cstdlib>

namespace {
    const int TABLE_SIZE = 512; // must be power of two
    float random_numbers[TABLE_SIZE];
    /* float random_numbers_2[TABLE_SIZE][TABLE_SIZE]; */

    const float PI = 4.0f * std::atan(1.0f);
    
    float CosInterp(float x, float below, float above)
    {
        const float f = 0.5f * (1.0f - std::cos(PI * x));
        return below + f * (above - below);
    }

    float DCosInterpDX(float x, float below, float above)
    {
        const float dfdx = 0.5f * PI * std::sin(PI * x);
        return dfdx * (above - below);
    }
    
    float Octave(float lambda, float x)
    {
        x /= lambda;

        const int below = int(std::floor(x));
        const int above = int(std::ceil(x));

        const float rbelow = random_numbers[below & (TABLE_SIZE - 1)];
        const float rabove = random_numbers[above & (TABLE_SIZE - 1)];

        return CosInterp(x - below, rbelow, rabove);
    }

    float DOctaveDX(float lambda, float x)
    {
        x /= lambda;

        const int below = int(std::floor(x));
        const int above = int(std::ceil(x));

        const float rbelow = random_numbers[below & (TABLE_SIZE - 1)];
        const float rabove = random_numbers[above & (TABLE_SIZE - 1)];

        return DCosInterpDX(x - below, rbelow, rabove) / lambda;
    }

    /*
    float Octave2(float lambda_x, float lambda_y, float x, float y)
    {
        x /= lambda_x;
        y /= lambda_y;

        const int below_x = int(std::floor(x));
        const int above_x = int(std::ceil(x));
        const int below_y = int(std::floor(y));
        const int above_y = int(std::ceil(y));

        const float rbxby = random_numbers_2[below_x & (TABLE_SIZE - 1)][below_y & (TABLE_SIZE - 1)];
        const float rbxay = random_numbers_2[below_x & (TABLE_SIZE - 1)][above_y & (TABLE_SIZE - 1)];
        const float raxby = random_numbers_2[above_x & (TABLE_SIZE - 1)][below_y & (TABLE_SIZE - 1)];
        const float raxay = random_numbers_2[above_x & (TABLE_SIZE - 1)][above_y & (TABLE_SIZE - 1)];

        const float rby = CosInterp(x - below_x, rbxby, raxby);
        const float ray = CosInterp(x - below_x, rbxay, raxay);
        
        return CosInterp(y - below_y, rby, ray);
    }
    */
}

void InitPerlin()
{
    static bool initialized = false;
    if (initialized) return;
    initialized = true;
    
    for (int i = 0; i < TABLE_SIZE; ++i) {
        random_numbers[i] = 2 * float(std::rand()) / float(RAND_MAX) - 1;
    }

    /*
    for (int i = 0; i < TABLE_SIZE; ++i) {
        for (int j = 0; j < TABLE_SIZE; ++j) {
            random_numbers_2[i][j] = 2 * float(std::rand()) / float(RAND_MAX) - 1;
        }
    }
    */
}

float Perlin(float lambda, float persistence, int octaves, float x)
{
    float amp = 1;
    float tot = 0;
    for (int i = 0; i < octaves; ++i) {
        tot += amp * Octave(lambda, x);
        lambda /= 2;
        x += 1997 * lambda;
        amp *= persistence;
    }
    return tot;
}

float DPerlinDX(float lambda, float persistence, int octaves, float x)
{
    float amp = 1;
    float tot = 0;
    for (int i = 0; i < octaves; ++i) {
        tot += amp * DOctaveDX(lambda, x);
        lambda /= 2;
        x += 1997 * lambda;
        amp *= persistence;
    }
    return tot;
}

/*
float Perlin2(float lambda_x, float lambda_y, float persistence, int octaves, float x, float y)
{
    float amp = 1;
    float tot = 0;
    for (int i = 0; i < octaves; ++i) {
        tot += amp * Octave2(lambda_x, lambda_y, x, y);
        lambda_x /= 2;
        lambda_y /= 2;
        x += 1997 * lambda_x;
        y += 1979 * lambda_y;
        amp *= persistence;
    }
    return tot;
}
*/
