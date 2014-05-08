/*
 * FILE:
 *   perlin.hpp
 *
 * PURPOSE:
 *   Simple Perlin noise
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

#ifndef PERLIN_HPP
#define PERLIN_HPP

void InitPerlin();

float Perlin(float base_lambda, float persistence, int octaves, float x);
float DPerlinDX(float base_lambda, float persistence, int octaves, float x);

/* 2D perlin noise, not used
float Perlin2(float base_lambda_x, float base_lambda_y, float persistence, int octaves, float x, float y);
*/

#endif
