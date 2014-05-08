/*
 * FILE:
 *   terrain_heightfield.hpp
 *
 * PURPOSE:
 *   Function to compute the terrain heightfield given current
 *   settings
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * CREATED:
 *   24-Oct-2011
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

#ifndef TERRAIN_HEIGHTFIELD_HPP
#define TERRAIN_HEIGHTFIELD_HPP

#include "boost/scoped_array.hpp"

void UpdateTerrainHeightfield();
float GetTerrainHeight(float x, float y);  // does interpolation / clamping

struct TerrainEntry {
    float B, dBdx, dBdy;
};

struct BottomEntry {
    float BY, BX, BA;
};

// used to initialize the heightfield texture (used by the terrain vertex shader).
// includes ghost zones.
// TODO: this should probably be baked into the terrain mesh instead.
extern boost::scoped_array<TerrainEntry> g_terrain_heightfield;

// used to initialize the "bottom" texture (used for simulation).
extern boost::scoped_array<BottomEntry> g_bottom;


extern float g_inlet_x;

#endif
