/*
 * FILE:
 *   terrain_heightfield.cpp
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

#include "perlin.hpp"
#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include <cmath>

namespace {
    float CalcV(float H, float W, float shape, float x)
    {
        return H * std::pow(fabs(x) / (W/2), shape);
    }

    float DVDX(float H, float W, float shape, float x)
    {
        float result = H * shape * std::pow(std::fabs(x) / (W/2), shape-1) / (W/2);
        if (x<0) result = -result;
        return result;
    }

    int RoundDown(float x)
    {
        return int(std::floor(x));
    }

    int RoundUp(float x)
    {
        return int(std::ceil(x));
    }

    float Lerp(float x, float down, float up)
    {
        return down + x * (up - down);
    }


    // height calculation.

    // for efficiency we split this API into three separate calls:
    // InitHeight() called once at the beginning
    // InitRow(float y) called once per y-value
    // GetHeight(float x, float &B, float &dB_dx, float &dB_dy) called once per point being sampled.

    float L, W, H, m_top, m_bottom, shape;
    float C_D_top, C_D_bottom, C_W_top, C_W_bottom;
    bool dam_on;
    float D_H, D_Y, D_MW, D_MH, D_T;
    float M_lambda, M_A, M_P;

    const int NUM_OCTAVES = 8;

    void InitHeight()
    {
        InitPerlin();

        L = GetSetting("valley_length");
        W = GetSetting("valley_width");
        H = GetSetting("valley_wall_height");
        m_top = GetSetting("gradient_top");
        m_bottom = GetSetting("gradient_bottom");
        shape = GetSetting("valley_shape");
        
        C_D_top = GetSetting("channel_depth_top");
        C_D_bottom = GetSetting("channel_depth_bottom");
        C_W_top = GetSetting("channel_width_top");
        C_W_bottom = GetSetting("channel_width_bottom");

        dam_on = GetIntSetting("dam_on") != 0;
        D_H = GetSetting("dam_height");
        D_Y = GetSetting("dam_position");
        D_MW = GetSetting("dam_middle_width");
        D_MH = GetSetting("dam_middle_height");
        D_T = GetSetting("dam_thickness");
        
        M_lambda = GetSetting("meander_wavelength");
        M_A = GetSetting("meander_amplitude");
        M_P = GetSetting("meander_fractal");
    }

    float z0, dz0_dy;
    float C_W, C_D, dCW_dy, dCD_dy;
    float R_C, dRC_dy, R_L, R_R, dRL_dy, dRR_dy;
    float B_L, B_R, dBL_dy, dBR_dy;
    float B_base, dBbase_dy;
    float dam_drop, ddamdrop_dy;
    
    void InitRow(float y)
    {
        // flatten outside of domain
        //y = std::max(0.0f, std::min(L, y));

        // reflect outside of domain
        const float dy = L / (GetSetting("mesh_size_y")-1);
        if (y < -dy/2) y = -dy - y;
        if (y > L+dy/2) y = 2*L + dy - y;


        z0 = y * (m_bottom + (m_top - m_bottom) * y / (2*L));
        dz0_dy = m_bottom + (m_top - m_bottom) * y / L;

        C_W = C_W_bottom + y * (C_W_top - C_W_bottom) / L;
        C_D = C_D_bottom + y * (C_D_top - C_D_bottom) / L;

        dCW_dy = (C_W_top - C_W_bottom) / L;
        dCD_dy = (C_D_top - C_D_bottom) / L;
        
        R_C = M_A * Perlin(M_lambda, M_P, NUM_OCTAVES, y);
        dRC_dy = M_A * DPerlinDX(M_lambda, M_P, NUM_OCTAVES, y);
        
        R_L = R_C - C_W/2;
        R_R = R_C + C_W/2;

        dRL_dy = dRC_dy - dCW_dy/2;
        dRR_dy = dRC_dy + dCW_dy/2;

        B_L = z0 + CalcV(H, W, shape, R_L);
        B_R = z0 + CalcV(H, W, shape, R_R);

        dBL_dy = dz0_dy + DVDX(H, W, shape, R_L) * dRL_dy;
        dBR_dy = dz0_dy + DVDX(H, W, shape, R_R) * dRR_dy;
        
        // B_base = min(B_L, B_R) - C_D
        if (B_L < B_R) {
            B_base = B_L - C_D;
            dBbase_dy = dBL_dy - dCD_dy;
        } else {
            B_base = B_R - C_D;
            dBbase_dy = dBR_dy - dCD_dy;
        }

        if (dam_on) {
            if (y < D_Y - D_T/2) {
                dam_drop = 2 * (D_Y - D_T/2 - y);
                ddamdrop_dy = -2;
            } else if (y > D_Y + D_T/2) {
                dam_drop = 2 * (y - D_Y - D_T/2);
                ddamdrop_dy = 2;
            } else {
                dam_drop = 0;
                ddamdrop_dy = 0;
            }
        } else {
            dam_drop = ddamdrop_dy = 0;
        }
    }

    void GetHeight(float x, float &B)
    {
        // reflect outside of domain
        const float dx = W / (GetSetting("mesh_size_x")-1);
        if (x < -W/2 - dx/2) x = -W - dx - x;
        if (x > W/2 + dx/2) x = W + dx - x;

        float B_star;

        if (x < R_L || x > R_R) {
            B_star = z0 + CalcV(H, W, shape, x);
            
        } else if (x < R_L + C_W/10) {
            B_star = B_L + (x - R_L) / (C_W/10) * (B_base - B_L);
            
        } else if (x < R_R - C_W/10) {
            B_star = B_base;
            
        } else {
            B_star = B_R + (R_R - x) / (C_W/10) * (B_base - B_R);
        }
        
        float B_D = B_base + D_H;
        const float T = std::fabs(D_MH)/2;
        const float S = D_MH < 0 ? -2.0f : 2.0f;
        
        if (x < R_C - D_MW/2 - T) {
            // do nothing
            
        } else if (x < R_C - D_MW/2) {
            B_D += S * (x - (R_C - D_MW/2 - T));
            
        } else if (x < R_C + D_MW/2) {
            B_D += D_MH;

        } else if (x < R_C + D_MW/2 + T) {
            B_D += S * ((R_C + D_MW/2 + T) - x);
        }

        B_D -= dam_drop;

        if (!dam_on || B_star > B_D) {
            B = B_star;
        } else {
            B = B_D;
        }
    }

    void GetHeightDeriv(float x, float &dB_dx, float &dB_dy)
    {
        // reflect outside of domain
        const float dx = W / (GetSetting("mesh_size_x")-1);
        if (x < -W/2 - dx/2) x = -W - dx - x;
        if (x > W/2 + dx/2) x = W + dx - x;

        float B, B_star, dBstar_dx, dBstar_dy;

        if (x < R_L || x > R_R) {
            B_star = z0 + CalcV(H, W, shape, x);
            dBstar_dx = DVDX(H, W, shape, x);
            dBstar_dy = dz0_dy;
            
        } else if (x < R_L + C_W/10) {
            B_star = B_L + (x - R_L) / (C_W/10) * (B_base - B_L);
            dBstar_dx = (B_base - B_L) / (C_W/10);
            
            const float deriv = -dRL_dy * (B_base - B_L) + (x - R_L) * (dBbase_dy - dBL_dy);
            dBstar_dy = dBL_dy + 10 * (C_W * deriv - (x - R_L) * (B_base - B_L) * dCW_dy) / (C_W * C_W);
            
        } else if (x < R_R - C_W/10) {
            B_star = B_base;
            dBstar_dx = 0;
            dBstar_dy = dBbase_dy;
            
        } else {
            B_star = B_R + (R_R - x) / (C_W/10) * (B_base - B_R);
            dBstar_dx = (B_R - B_base) / (C_W/10);
            
            const float deriv = dRR_dy * (B_base - B_R) + (R_R - x) * (dBbase_dy - dBR_dy);
            dBstar_dy = dBR_dy + 10 * (C_W * deriv - (R_R - x) * (B_base - B_R) * dCW_dy) / (C_W * C_W);
        }
        
        float B_D = B_base + D_H;
        float dBD_dx = 0;
        const float T = std::fabs(D_MH)/2;
        const float S = D_MH < 0 ? -2.0f : 2.0f;
        
        if (x < R_C - D_MW/2 - T) {
            // do nothing
            
        } else if (x < R_C - D_MW/2) {
            B_D += S * (x - (R_C - D_MW/2 - T));
            dBD_dx = S;
            
        } else if (x < R_C + D_MW/2) {
            B_D += D_MH;

        } else if (x < R_C + D_MW/2 + T) {
            B_D += S * ((R_C + D_MW/2 + T) - x);
            dBD_dx = -S;
        }

        B_D -= dam_drop;

        if (!dam_on || B_star > B_D) {
            B = B_star;
            dB_dx = dBstar_dx;
            dB_dy = dBstar_dy;
        } else {
            B = B_D;
            dB_dx = dBD_dx;
            dB_dy = dBbase_dy - ddamdrop_dy;
        }
    }

}

boost::scoped_array<TerrainEntry> g_terrain_heightfield;
boost::scoped_array<BottomEntry> g_bottom;
float g_inlet_x;

float VH = 10.0;
void GetHeightDeriv2(float x, float y, float&dBdx, float& dBdy)
{
    float fx = std::fabs(x);
    dBdy = 0.0;
    //dBdx = ( (x / (W/2)) * VH ) / (fx);
    dBdx = 0.0;//( (x / (W/2)) * VH ) / (fx);

}

void GetHeight2(float x, float y, float& h)
{
    if( x < 0.0 && y < 100.0 )
    {
        float fx = std::fabs(x);
        h = (fx / (W/2)) * VH;// + 10.0f;

    }
    else
    {
        float fx = std::fabs(x);
        h = (fx / (W/2)) * VH;
    }
}

void UpdateTerrainHeightfield()
{
    using namespace std;

    const int nx = GetIntSetting("mesh_size_x");//300
    const int ny = GetIntSetting("mesh_size_y");//900
    const int pitch = nx+4;
    
    g_terrain_heightfield.reset(new TerrainEntry[pitch * (ny+4)]);
    g_bottom.reset(new BottomEntry[pitch * (ny+4)]);
    
    //InitHeight();
    L = GetSetting("valley_length");//100 -50 -> 50
    W = GetSetting("valley_width");//300 0 -> 300


    // Calculate B at each mesh point (cell centre)
    for (int j = 0; j < ny + 4; ++j) {
        //const float y = float(j-2) / float(ny-1) * L;
        const float y = L - float(j-2) / float(ny-1) * L;

        //InitRow(y);
        
        for (int i = 0; i < nx + 4; ++i) {
            //const float x = float(i-2) / float(nx-1) * W - (W/2);
            const float x = (W)*0.5 - float(i-2) / float(nx-1) * W;
            
            TerrainEntry &out = g_terrain_heightfield[j * pitch + i];
            //GetHeightDeriv(x, out.dBdx, out.dBdy);
            GetHeightDeriv2(x, y, out.dBdx, out.dBdy);
            out.dBdx = out.dBdx;
        }
    }

    
    // BX(i,j) = 0.5 * (B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    // BY(i,j) = 0.5 * (B(i-1/2, j+1/2) + B(i+1/2, j+1/2))
    // BA(i,j) = 0.25 * (B(i-1/2, j-1/2) + B(i-1/2, j+1/2) + B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    
    for (int j = 0; j < ny + 5; ++j) {
        const float y = L - (float(j)-2.5f) / float(ny-1) * L;    // j - 1/2
        //const float y = (float(j)-2.5f) / float(ny-1) * L;    // j - 1/2
        //InitRow(y);

        for (int i = 0; i < nx + 5; ++i) {
            //const float x = (float(i)-2.5f) / float(nx-1) * W - (W/2);   // i - 1/2
            const float x = (W)*0.5 - (float(i)-2.5f) / float(nx-1) * W;   // i - 1/2

            float h;
            //GetHeight(x, h);
            GetHeight2(x, y, h);

            if (i < nx+4) {   // (i)
                if (j < ny+4) {   // (i,j)
                    g_bottom[j*pitch + i].BA = 0.25f * h;
                }
                if (j > 0) {   // (i,j-1)
                    g_bottom[(j-1)*pitch + i].BA += 0.25f * h;
                    g_bottom[(j-1)*pitch + i].BY = 0.5f * h;
                }
            }

            if (i > 0) {  // (i-1)
                if (j < ny+4) {   // (i-1,j)
                    g_bottom[j*pitch + i-1].BA += 0.25f * h;
                    g_bottom[j*pitch + i-1].BX = 0.5f * h;
                }

                if (j > 0) {  // (i-1,j-1)
                    g_bottom[(j-1)*pitch + i-1].BA += 0.25f * h;
                    g_bottom[(j-1)*pitch + i-1].BX += 0.5f * h;
                    g_bottom[(j-1)*pitch + i-1].BY += 0.5f * h;
                }
            }
        }
    }

    // We now use BA instead of B in g_terrain_heightfield.
    // This prevents water "showing through" in steep areas.
    for (int j = 0; j < ny + 4; ++j) {
        for (int i = 0; i < nx + 4; ++i) {
            g_terrain_heightfield[j * pitch + i].B = g_bottom[j * pitch + i].BA;
        }
    }

    g_inlet_x = GetSetting("meander_amplitude") * 
        Perlin(GetSetting("meander_wavelength"), GetSetting("meander_fractal"), NUM_OCTAVES, L);
}

void UpdateTerrainHeightfield1()
{
    using namespace std;

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const int pitch = nx+4;
    
    g_terrain_heightfield.reset(new TerrainEntry[pitch * (ny+4)]);
    g_bottom.reset(new BottomEntry[pitch * (ny+4)]);
    
    InitHeight();

    // Calculate B at each mesh point (cell centre)
    for (int j = 0; j < ny + 4; ++j) {
        const float y = float(j-2) / float(ny-1) * L;

        InitRow(y);
        
        for (int i = 0; i < nx + 4; ++i) {
            const float x = float(i-2) / float(nx-1) * W - (W/2);
            
            TerrainEntry &out = g_terrain_heightfield[j * pitch + i];
            GetHeightDeriv(x, out.dBdx, out.dBdy);
        }
    }

    
    // BX(i,j) = 0.5 * (B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    // BY(i,j) = 0.5 * (B(i-1/2, j+1/2) + B(i+1/2, j+1/2))
    // BA(i,j) = 0.25 * (B(i-1/2, j-1/2) + B(i-1/2, j+1/2) + B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    
    for (int j = 0; j < ny + 5; ++j) {
        const float y = (float(j)-2.5f) / float(ny-1) * L;    // j - 1/2
        InitRow(y);

        for (int i = 0; i < nx + 5; ++i) {
            const float x = (float(i)-2.5f) / float(nx-1) * W - (W/2);   // i - 1/2

            float h;
            GetHeight(x, h);

            if (i < nx+4) {   // (i)
                if (j < ny+4) {   // (i,j)
                    g_bottom[j*pitch + i].BA = 0.25f * h;
                }
                if (j > 0) {   // (i,j-1)
                    g_bottom[(j-1)*pitch + i].BA += 0.25f * h;
                    g_bottom[(j-1)*pitch + i].BY = 0.5f * h;
                }
            }

            if (i > 0) {  // (i-1)
                if (j < ny+4) {   // (i-1,j)
                    g_bottom[j*pitch + i-1].BA += 0.25f * h;
                    g_bottom[j*pitch + i-1].BX = 0.5f * h;
                }

                if (j > 0) {  // (i-1,j-1)
                    g_bottom[(j-1)*pitch + i-1].BA += 0.25f * h;
                    g_bottom[(j-1)*pitch + i-1].BX += 0.5f * h;
                    g_bottom[(j-1)*pitch + i-1].BY += 0.5f * h;
                }
            }
        }
    }

    // We now use BA instead of B in g_terrain_heightfield.
    // This prevents water "showing through" in steep areas.
    for (int j = 0; j < ny + 4; ++j) {
        for (int i = 0; i < nx + 4; ++i) {
            g_terrain_heightfield[j * pitch + i].B = g_bottom[j * pitch + i].BA;
        }
    }

    g_inlet_x = GetSetting("meander_amplitude") * 
        Perlin(GetSetting("meander_wavelength"), GetSetting("meander_fractal"), NUM_OCTAVES, L);
}

void UpdateTerrainHeightfield2()
{
    using namespace std;

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const int pitch = nx+4;
    
    g_terrain_heightfield.reset(new TerrainEntry[pitch * (ny+4)]);
    g_bottom.reset(new BottomEntry[pitch * (ny+4)]);
    
    InitHeight();

    // Calculate B at each mesh point (cell centre)
    for (int j = 0; j < ny + 4; ++j) {
        const float y = float(j-2) / float(ny-1) * L;

        InitRow(y);
        
        for (int i = 0; i < nx + 4; ++i) {
            const float x = float(i-2) / float(nx-1) * W - (W/2);
            
            TerrainEntry &out = g_terrain_heightfield[j * pitch + i];
            GetHeightDeriv(x, out.dBdx, out.dBdy);
            //out.dBdx = 1.0;
            //out.dBdy = 1.0;
            g_bottom[j * pitch + i].BA = 0.0;
            g_bottom[j * pitch + i].BX = 0.0;
            g_bottom[j * pitch + i].BY = 0.0;
            
        }
    }

    
    // BX(i,j) = 0.5 * (B(i+1/2, j-1/2) + B(i+1/2, j+1/2))
    // BY(i,j) = 0.5 * (B(i-1/2, j+1/2) + B(i+1/2, j+1/2))
    // BA(i,j) = 0.25 * (B(i-1/2, j-1/2) + B(i-1/2, j+1/2) + B(i+1/2, j-1/2) + B(i+1/2, j+1/2))

    // We now use BA instead of B in g_terrain_heightfield.
    // This prevents water "showing through" in steep areas.
    for (int j = 0; j < ny + 4; ++j) {
        for (int i = 0; i < nx + 4; ++i) {
            g_terrain_heightfield[j * pitch + i].B = g_bottom[j * pitch + i].BA;
        }
    }

    g_inlet_x = GetSetting("meander_amplitude") * 
        Perlin(GetSetting("meander_wavelength"), GetSetting("meander_fractal"), NUM_OCTAVES, L);
}

float GetTerrainHeight(float x, float y)
{
    const int width = GetIntSetting("mesh_size_x");
    const int height = GetIntSetting("mesh_size_y");    
    
    const float L = GetSetting("valley_length");
    const float W = GetSetting("valley_width");

    // convert x,y to [0,width-1], [0,height-1] scale
    x = (x + W/2) / W * (width-1);
    y = y / L * (height-1);

    if (x < 0) x = 0;
    if (x > width-1) x = float(width-1);
    if (y < 0) y = 0;
    if (y > height-1) y = float(height-1);

    // now add 2 to account for ghost cells
    x += 2;
    y += 2;
    
    const int xdown = RoundDown(x);
    const int xup = RoundUp(x);
    const int ydown = RoundDown(y);
    const int yup = RoundUp(y);

    const float xfrac = x - xdown;
    const float yfrac = y - ydown;

    const int pitch = width + 4;
    
    const float result_y_down = Lerp(xfrac,
                                     g_terrain_heightfield[ydown * pitch + xdown].B,
                                     g_terrain_heightfield[ydown * pitch + xup].B);
    const float result_y_up = Lerp(xfrac,
                                   g_terrain_heightfield[yup * pitch + xdown].B,
                                   g_terrain_heightfield[yup * pitch + xup].B);

    return Lerp(yfrac, result_y_down, result_y_up);
}
