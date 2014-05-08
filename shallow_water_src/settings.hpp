/*
 * FILE:
 *   settings.hpp
 *
 * PURPOSE:
 *   Abstract interface for rendering and simulation settings.
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

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <vector>

enum SettingType {
    S_NULL,
    S_NEW_TAB,
    S_LABEL,
    S_SLIDER,
    S_SLIDER_INT,
    S_SLIDER_MULT_4,
    S_CHECKBOX,
    S_LEFT_MOUSE_DROPDOWN
};

enum ResetType {
    R_NONE,
    R_TERRAIN,
    R_MESH,     // remesh, no water
    R_VALLEY,   // remesh, setup water in valley
    R_SEA,      // remesh, fill to sea level
    R_SQUARE    // remesh, fill in a square in the middle
};

enum LeftMouseSettings {
    LM_ADD_WATER,
    LM_REMOVE_WATER,
    LM_STIR_WATER,
    LM_RAISE_TERRAIN,
    LM_LOWER_TERRAIN
};

struct Setting {
    const char * name;
    const char * unit;
    SettingType type;
    ResetType reset_type;
    double min;
    double max;
    double value;
};

// global array of settings, 'null terminated'
extern Setting g_settings[];

// flag to communicate when the settings have changed
extern ResetType g_reset_type;

// helper functions
float GetSetting(const char * name);
inline float GetSetting(const std::string &s) { return GetSetting(s.c_str()); }
int GetIntSetting(const char *name);

void SetSetting(const char *name, float new_value);
inline void SetSettingD(const char *name, double val) { SetSetting(name, float(val)); }

#endif
