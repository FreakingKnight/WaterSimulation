/*
 * FILE:
 *   settings.cpp
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

#include "settings.hpp"

#include <cstdlib>
#include <cmath>
#include <map>
#include <string>

namespace {
    const double PI = std::atan(1.0) * 4.0;

    class CompareStr {
    public:
        bool operator()(const char *lhs, const char *rhs) const {
            return strcmp(lhs, rhs) < 0;
        }
    };
    
    std::map<const char*, Setting*, CompareStr> g_name_to_setting;
}

Setting g_settings[] =
    {
        // name, unit, type, min, max, default


        // SETTINGS TAB
        
        { "settings", "", S_NEW_TAB },
        
        { "mesh_size_x", "", S_SLIDER_MULT_4, R_MESH, 10, 1200 },
        { "mesh_size_y", "", S_SLIDER_MULT_4, R_MESH, 10, 1200 },
        {""},
        
        { "solid_walls", "", S_CHECKBOX, R_NONE, 0, 1 },
        { "inflow_width", "m", S_SLIDER, R_NONE, 0, 100 },
        { "inflow_height", "m", S_SLIDER, R_NONE, 0, 10 },
        { "" },

        { "gravity", "m s^-2", S_SLIDER, R_NONE, 0, 30 },
        { "friction", "m s^-1", S_SLIDER, R_NONE, 0, 0.5 },
        { "" },

        { "theta", "", S_SLIDER, R_NONE, 1, 2 },
        { "max_cfl_number", "", S_SLIDER, R_NONE, 0, 0.6 },
        { "timesteps_per_frame", "", S_SLIDER_INT, R_NONE, 1, 10 },
        { "time_acceleration", "", S_SLIDER, R_NONE, 1, 30 },
        {""},

        { "clip_camera", "", S_CHECKBOX, R_NONE, 0, 1 },



        // TERRAIN TAB

        { "terrain", "", S_NEW_TAB },

        { "valley_length", "m", S_SLIDER, R_TERRAIN, 10, 1000 },
        { "valley_width", "m", S_SLIDER, R_TERRAIN, 10, 1000 },
        { "valley_wall_height", "m", S_SLIDER, R_TERRAIN, 0, 100 },
        { "gradient_top", "", S_SLIDER, R_TERRAIN, 0, 0.5 },
        { "gradient_bottom", "", S_SLIDER, R_TERRAIN, 0, 0.5 },
        { "valley_shape", "", S_SLIDER, R_TERRAIN, 1, 4 },
        { "" },
        
        { "channel_depth_top", "m", S_SLIDER, R_TERRAIN, 0, 10 },
        { "channel_depth_bottom", "m", S_SLIDER, R_TERRAIN, 0, 10 },
        { "channel_width_top", "m", S_SLIDER, R_TERRAIN, 0, 100 },
        { "channel_width_bottom", "m", S_SLIDER, R_TERRAIN, 0, 100 },
        { "" },
        
        { "dam_on", "", S_CHECKBOX, R_TERRAIN, 0, 1 },
        { "dam_height", "m", S_SLIDER, R_TERRAIN, 0, 100 },
        { "dam_position", "m", S_SLIDER, R_TERRAIN, 0, 1000 },
        { "dam_middle_width", "m", S_SLIDER, R_TERRAIN, 0, 50 },
        { "dam_middle_height", "m", S_SLIDER, R_TERRAIN, -10, 10 },
        { "dam_thickness", "m", S_SLIDER, R_TERRAIN, 0, 100 },
        { "" },
        
        { "meander_wavelength", "m", S_SLIDER, R_TERRAIN, 10, 1000 },
        { "meander_amplitude", "m", S_SLIDER, R_TERRAIN, 0, 50 },
        { "meander_fractal", "", S_SLIDER, R_TERRAIN, 0, 1 },


        // SEA TAB
        
        { "sea", "", S_NEW_TAB },

        { "use_sea_level", "", S_CHECKBOX, R_NONE, 0, 1 },
        { "sea_level", "m", S_SLIDER, R_NONE, 1, 25 },
        { "" },
        
        { "sa1", "m",    S_SLIDER, R_NONE, 0, 1 },
        { "sk1", "m^-1", S_SLIDER, R_NONE, 0, 1 },
        { "sk1_dir", "", S_SLIDER, R_NONE, 0, 2*PI },
        { "so1", "s^-1", S_SLIDER, R_NONE, 0, 10 },
        { "" },
        
        { "sa2", "m",    S_SLIDER, R_NONE, 0, 1 },
        { "sk2", "m^-1", S_SLIDER, R_NONE, 0, 1 },
        { "sk2_dir", "", S_SLIDER, R_NONE, 0, 2*PI },
        { "so2", "s^-1", S_SLIDER, R_NONE, 0, 10 },
        { "" },
        
        { "sa3", "m",    S_SLIDER, R_NONE, 0, 1 },
        { "sk3", "m^-1", S_SLIDER, R_NONE, 0, 1 },
        { "sk3_dir", "", S_SLIDER, R_NONE, 0, 2*PI },
        { "so3", "s^-1", S_SLIDER, R_NONE, 0, 10 },
        { "" },
        
        { "sa4", "m",    S_SLIDER, R_NONE, 0, 1 },
        { "sk4", "m^-1", S_SLIDER, R_NONE, 0, 1 },
        { "sk4_dir", "", S_SLIDER, R_NONE, 0, 2*PI },
        { "so4", "s^-1", S_SLIDER, R_NONE, 0, 10 },
        { "" },


        // GRAPHICS TAB

        { "graphics", "", S_NEW_TAB },

        { "fov", "deg", S_SLIDER, R_NONE, 1, 150 },
        {""},
        { "sun_alt", "deg", S_SLIDER, R_NONE, 0, 90 },
        { "sun_az", "deg", S_SLIDER, R_NONE, 0, 360 },
        { "ambient", "", S_SLIDER, R_NONE, 0, 2 },
        { "" },
        { "fresnel_coeff", "", S_SLIDER, R_NONE, 0, 1 },
        { "fresnel_exponent", "", S_SLIDER, R_NONE, 0, 10 },
        { "specular_intensity", "", S_SLIDER, R_NONE, 0, 2 },
        { "specular_exponent", "", S_SLIDER, R_NONE, 0, 60 },
        { "refractive_index", "", S_SLIDER, R_NONE, 1, 3 },
        { "attenuation_1", "", S_SLIDER, R_NONE, 0, 0.5 },
        { "attenuation_2", "", S_SLIDER, R_NONE, 0, 2.0 },
        { "deep_r", "", S_SLIDER, R_NONE, 0, 1 },
        { "deep_g", "", S_SLIDER, R_NONE, 0, 1 },
        { "deep_b", "", S_SLIDER, R_NONE, 0, 1 },


        // NON TABBED WIDGETS

        { "", "", S_NEW_TAB },

        { "left_mouse_action", "", S_LEFT_MOUSE_DROPDOWN, R_NONE, 0, 4, 0 },
        { "left_mouse_radius", "", S_SLIDER, R_NONE, 1, 25, 5 },
        { "left_mouse_strength", "", S_SLIDER, R_NONE, 1, 50, 10 },
        { "" },
        
        { "mass", "kg", S_LABEL },
        { "x_momentum", "kg m s^-1", S_LABEL },
        { "y_momentum", "kg m s^-1", S_LABEL },
        { "kinetic_energy", "J", S_LABEL },
        { "potential_energy", "J", S_LABEL },
        { "total_energy", "J", S_LABEL },
        { "testable", "J", S_LABEL },
        { "testable1", "J", S_LABEL },
        { "max_speed", "m s^-1", S_LABEL },
        { "max_depth", "m", S_LABEL },
        { "max_froude_number", "", S_LABEL },
        { "" },

        { "fps", "s^-1", S_LABEL },
        { "timestep", "s", S_LABEL },
        { "cfl_number", "", S_LABEL },
        { "time_ratio", "", S_LABEL },

        // mouse position information
        {""},
        { "x", "", S_LABEL },
        { "y", "", S_LABEL },
        { "z", "", S_LABEL },
        { "depth", "", S_LABEL },
        { "u", "", S_LABEL },
        { "v", "", S_LABEL },
        { "froude", "", S_LABEL },
        
        {}   // Null terminator
    };

ResetType g_reset_type;

namespace {
    void InitNameToSetting()
    {
        for (Setting *p = &g_settings[0]; p->name; ++p) {
            g_name_to_setting.insert(std::make_pair(p->name, p));
        }
    }
}

float GetSetting(const char *name)
{
    if (g_name_to_setting.empty()) InitNameToSetting();

    std::map<const char*, Setting*, CompareStr>::iterator it = g_name_to_setting.find(name);
    if (it == g_name_to_setting.end()) std::abort();
    
    return float(it->second->value);
}

int GetIntSetting(const char *name)
{
    return int(GetSetting(name) + 0.5f);
}

void SetSetting(const char *name, float new_value)
{
    if (g_name_to_setting.empty()) InitNameToSetting();
    std::map<const char*, Setting*, CompareStr>::iterator it = g_name_to_setting.find(name);
    if (it != g_name_to_setting.end()) it->second->value = double(new_value);
}
