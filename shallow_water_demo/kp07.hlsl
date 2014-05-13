/* -*- c++ -*-
 * 
 * FILE:
 *   kp07.hlsl
 *
 * PURPOSE:
 *   HLSL implementation of the Kurganov-Petrova scheme for solving
 *   the 2D shallow water equations. See:
 *
 *   A. Kurganov and G. Petrova, "A Second-Order Well-Balanced
 *   Positivity Preserving Central-Upwind Scheme for the Saint-Venant
 *   System", Commun. Math. Sci. 5, 133-160, 2007.
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * CREATED:
 *   27-Oct-2011
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

cbuffer SimConstBuffer : register( b0 )
{
    // THETA = parameter for minmod limiter. between 1 and 2 inclusive.
    // 1 = more dissipative, 2 = more oscillatory. 1.3 is a good default.
    // (Note: we actually store 2*THETA.)
    float TWO_THETA;

    float two_over_nx_plus_four;
    float two_over_ny_plus_four;

    float g;
    float half_g;
    float g_over_dx;
    float g_over_dy;
    float one_over_dx;
    float one_over_dy;
    float dt;
    float epsilon;      // suggestion: dx^4 or dy^4. but that may be too small for floats?

    int nx, ny;  // number of cells, excluding ghost zones

    float friction;              // m s^-1
};

cbuffer BoundaryConstBuffer : register( b0 )
{
    float boundary_epsilon;
    int reflect_x, reflect_y;
    int solid_wall_flag;
    int inflow_x_min, inflow_x_max;
    float sea_level, inflow_height, inflow_speed;
    float boundary_g;
    float total_time;
    float sa1, skx1, sky1, so1;
    float sa2, skx2, sky2, so2;
    float sa3, skx3, sky3, so3;
    float sa4, skx4, sky4, so4;
    float sdecay;
};

// .r = B(j, k+1/2)   "BN" or "BY"
// .g = B(j+1/2, k)   "BE" or "BX"
// .b = B(j, k)       "BA"
// Note: this is the same size as the main textures i.e. includes ghost zones.
Texture2D<float3> txBottom : register( t1 );

// .r = w_bar
// .g = hu_bar
// .b = hv_bar
// .a = (unused)
Texture2D<float4> txState : register( t0 );

// .r = N
// .g = E
// .b = S
// .a = W
Texture2D<float4> txH : register( t0 );
Texture2D<float4> txU : register( t1 );
Texture2D<float4> txV : register( t2 );

// .r = 1 (w-flux)
// .g = 2 (hu-flux)
// .b = 3 (hv-flux)
// .a = (unused)
Texture2D<float4> txXFlux : register( t2 );
Texture2D<float4> txYFlux : register( t3 );

struct VS_INPUT {
    float2 tex_idx : TEX_IDX;
};

struct VS_OUTPUT {
    float4 pos : SV_POSITION;
    float2 tex_idx : TEX_IDX; // should be cast to integer.
};

struct PASS_1_OUTPUT {
    float4 h : SV_TARGET0;    // {hN, hE, hS, hW}
    float4 u : SV_TARGET1;    // {uN, uE, uS, uW}
    float4 v : SV_TARGET2;    // {vN, vE, vS, vW}
    float4 n : SV_TARGET3;    // {nX, nY, nZ, unused}
};

struct PASS_2_OUTPUT {
    float4 xflux : SV_TARGET0;   // {Hx1, Hx2, Hx3, unused}
    float4 yflux : SV_TARGET1;   // {Hy1, Hy2, Hy3, unused}
};

VS_OUTPUT SimVertexShader(VS_INPUT input)
{
    VS_OUTPUT output;
    output.pos.x = input.tex_idx.x * two_over_nx_plus_four - 1;
    output.pos.y = 1 - input.tex_idx.y * two_over_ny_plus_four;
    output.pos.z = 0.5f;
    output.pos.w = 1.0f;
    output.tex_idx = input.tex_idx;
    return output;
}

int3 GetTexIdx(VS_OUTPUT input)
{
    return int3(input.tex_idx.x, input.tex_idx.y, 0);
}

float MinMod(float a, float b, float c)
{
    return (a > 0 && b > 0 && c > 0) ? min(min(a,b),c)
        : (a < 0 && b < 0 && c < 0) ? max(max(a,b),c) : 0;
}

void Reconstruct(float west, float here, float east,
                 out float out_west, out float out_east)
{
    // west, here, east = values of U_bar at j-1, j, j+1 (or k-1, k, k+1)
    // out_west, out_east = reconstructed values of U_west and U_east at (j,k)
    
    float dx_grad_over_two = 0.25f * MinMod(TWO_THETA * (here - west),
                                            (east - west),
                                            TWO_THETA * (east - here));

    out_east = here + dx_grad_over_two;
    out_west = here - dx_grad_over_two;
}

void CorrectW(float B_west, float B_east, float w_bar,
              inout float w_west, inout float w_east)
{
    if (w_east < B_east) {
        w_east = B_east;
        w_west = max(B_west, 2 * w_bar - B_east);
            
    } else if (w_west < B_west) {
        w_east = max(B_east, 2 * w_bar - B_west);
        w_west = B_west;
    }
}              

void CalcUV(float4 h, float4 hu, float4 hv, out float4 u, out float4 v)
{
    // in:  {hN, hE, hS, hW},  {huN, huE, huS, huW},  {hvN, hvE, hvS, hvW}
    // out: {uN, uE, uS, uW},  {vN, vE, vS, vW}
    const float4 h2 = h * h;
    const float4 h4 = h2 * h2;
    const float4 divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, epsilon));
    u = divide_by_h * hu;
    v = divide_by_h * hv;
}

float CalcUV_Scalar(float h, float hu, float hv, out float u, out float v)
{
    const float h2 = h * h;
    const float h4 = h2 * h2;
    const float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, epsilon));
    u = divide_by_h * hu;
    v = divide_by_h * hv;
    return divide_by_h;
}

float CalcU_Boundary(float h, float hu)
{
    const float h2 = h * h;
    const float h4 = h2 * h2;
    const float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + max(h4, boundary_epsilon));
    return divide_by_h * hu;
}

float NumericalFlux(float aplus, float aminus, float Fplus, float Fminus, float Udifference)
{
    if (aplus - aminus > 0) {
        return (aplus * Fminus - aminus * Fplus + aplus * aminus * Udifference) / (aplus - aminus);
    } else {
        return 0;
    }
}


// Pass 1 -- Reconstruct h, u, v at the four edges of each cell.

// t0: txState
// t1: txBottom
// Output: txH, txU, txV, txNormal

// Runs on bulk + first ghost layer either side

PASS_1_OUTPUT Pass1(VS_OUTPUT input)
{
    // Read in relevant texture values

    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0);

    // in = {w, hu, hv, _} (cell average)
    const float4 in_here = txState.Load(idx);
    const float4 in_south = txState.Load(idx + int3(0, -1, 0));
    const float4 in_north = txState.Load(idx + int3(0, 1, 0));
    const float4 in_west = txState.Load(idx + int3(-1, 0, 0));
    const float4 in_east = txState.Load(idx + int3(1, 0, 0));

    float4 B;     // {BN, BE, BS, BW}
    B.rg = txBottom.Load(idx).rg;
    B.b = txBottom.Load(idx + int3(0, -1, 0)).r;
    B.a = txBottom.Load(idx + int3(-1, 0, 0)).g;
    
    
    // Reconstruct w, hu and hv at the four cell edges (N, E, S, W)
    
    float4 w;    // {wN, wE, wS, wW}
    float4 hu;   // {huN, huE, huS, huW}
    float4 hv;   // {hvN, hvE, hvS, hvW}
    
    Reconstruct(in_west.r, in_here.r, in_east.r, w.a, w.g);
    Reconstruct(in_south.r, in_here.r, in_north.r, w.b, w.r);
    
    Reconstruct(in_west.g, in_here.g, in_east.g, hu.a, hu.g);
    Reconstruct(in_south.g, in_here.g, in_north.g, hu.b, hu.r);
    
    Reconstruct(in_west.b, in_here.b, in_east.b, hv.a, hv.g);
    Reconstruct(in_south.b, in_here.b, in_north.b, hv.b, hv.r);


    // Correct the w values to ensure positivity of h

    CorrectW(B.a, B.g, in_here.r, w.a, w.g);    // wW and wE (from BW, BE, wbar)
    CorrectW(B.b, B.r, in_here.r, w.b, w.r);    // wS and wN (from BS, BN, wbar)


    // Reconstruct h from (corrected) w
    // Calculate u and v from h, hu and hv

    PASS_1_OUTPUT output;
    output.h = w - B;
    CalcUV(output.h, hu, hv, output.u, output.v);


    // Calculate normal 

    float3 normal;
    //normal.x = (w.g - w.a) * one_over_dx;
    //normal.y = (w.r - w.b) * one_over_dy;
    //normal.z = -1;
    normal.x = (in_west.r - in_east.r) * one_over_dx;
    normal.y = (in_south.r - in_north.r) * one_over_dy;
    normal.z = 2;
    normal = normalize(normal);
    output.n = float4(normal.x, normal.y, normal.z, 0);
    
    return output;
}


// Pass 2 -- Calculate fluxes

// t0: txH
// t1: txU
// t2: txV
// Output: txXFlux and txYFlux

// Runs on bulk + first ghost layer to west and south only

PASS_2_OUTPUT Pass2( VS_OUTPUT input )
{
    const int3 idx = GetTexIdx(input);
    
    const float2 h_here = txH.Load(idx).rg;                 // {hN, hE}   evaluated here
    const float hW_east = txH.Load(idx + int3(1,0,0)).a;    // hW evaluated at (j+1, k)
    const float hS_north = txH.Load(idx + int3(0,1,0)).b;   // hS evaluated at (j, k+1)

    const float2 u_here = txU.Load(idx).rg;
    const float uW_east = txU.Load(idx + int3(1,0,0)).a;
    const float uS_north = txU.Load(idx + int3(0,1,0)).b;

    const float2 v_here = txV.Load(idx).rg;
    const float vW_east = txV.Load(idx + int3(1,0,0)).a;
    const float vS_north = txV.Load(idx + int3(0,1,0)).b;
    
    // compute wave speeds
    const float2 cNE = sqrt(max(0, g * h_here.rg));    // {cN, cE} evaluated here
    const float cW = sqrt(max(0, g * hW_east));        // cW evaluated at (j+1, k)
    const float cS = sqrt(max(0, g * hS_north));       // cS evaluated at (j, k+1)

    // compute propagation speeds
    const float aplus  = max(max(u_here.g + cNE.g, uW_east + cW), 0);
    const float aminus = min(min(u_here.g - cNE.g, uW_east - cW), 0);
    const float bplus  = max(max(v_here.r + cNE.r, vS_north + cS), 0);
    const float bminus = min(min(v_here.r - cNE.r, vS_north - cS), 0);

    // compute fluxes
    PASS_2_OUTPUT output = (PASS_2_OUTPUT) 0;
    output.xflux.r = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * uW_east,
                                   h_here.g * u_here.g,
                                   hW_east - h_here.g);
    
    output.xflux.g = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * (uW_east * uW_east + half_g * hW_east),
                                   h_here.g * (u_here.g * u_here.g + half_g * h_here.g),
                                   hW_east * uW_east - h_here.g * u_here.g);

    output.xflux.b = NumericalFlux(aplus,
                                   aminus,
                                   hW_east * uW_east * vW_east,
                                   h_here.g * u_here.g * v_here.g,
                                   hW_east * vW_east - h_here.g * v_here.g);

    output.yflux.r = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * vS_north,
                                   h_here.r * v_here.r,
                                   hS_north - h_here.r);

    output.yflux.g = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * uS_north * vS_north,
                                   h_here.r * u_here.r * v_here.r,
                                   hS_north * uS_north - h_here.r * u_here.r);

    output.yflux.b = NumericalFlux(bplus,
                                   bminus,
                                   hS_north * (vS_north * vS_north + half_g * hS_north),
                                   h_here.r * (v_here.r * v_here.r + half_g * h_here.r),
                                   hS_north * vS_north - h_here.r * v_here.r);

    return output;
}


// Pass 3 -- Do timestep and calculate new w_bar, hu_bar, hv_bar.

// t0: txState
// t1: txBottom
// t2: txXFlux
// t3: txYFlux
// Output: New txState

// Runs on interior points only

float FrictionCalc(float h, float u)
{
    return max(h*u*dt*0.2f, friction * h * abs(u) * u);
}

float4 Pass3( VS_OUTPUT input ) : SV_Target
{
    const int3 idx = GetTexIdx(input);
    const float3 B_here = txBottom.Load(idx);

    float4 result;
    
    const float3 in_state = txState.Load(idx).rgb;   // w, hu and hv (cell avgs, evaluated here)
        
    const float3 xflux_here = txXFlux.Load(idx).rgb;
    const float3 xflux_west = txXFlux.Load(idx + int3(-1,0,0)).rgb;
    const float3 yflux_here = txYFlux.Load(idx).rgb;
    const float3 yflux_south = txYFlux.Load(idx + int3(0,-1,0)).rgb;
        
    const float BX_west = txBottom.Load(idx + int3(-1,0,0)).g;
    const float BY_south = txBottom.Load(idx + int3(0,-1,0)).r;

    // friction calculation
    const float h = max(0, in_state.r - B_here.b);
    float u, v;
    CalcUV_Scalar(h, in_state.g, in_state.b, u, v);

    const float3 source_term = 
        float3(0,
               -g_over_dx * h * (B_here.g - BX_west)   - FrictionCalc(h, u),
               -g_over_dy * h * (B_here.r - BY_south)  - FrictionCalc(h, v));
        
    const float3 d_by_dt =
        (xflux_west - xflux_here) * one_over_dx
        + (yflux_south - yflux_here) * one_over_dy
        + source_term;
        
    // simple Euler time stepping
    const float3 result3 = in_state + d_by_dt * dt;
    result = float4(result3.r, result3.g, result3.b, 0);

    return result;
}


// Boundary condition shaders

// fixed depth boundary calculation for east/north boundary
// returns h and hu for ghost zone (hv_ghost = 0).
void FixedHBoundary(float h_desired,
                    float h_real, float hu_real,
                    out float h_ghost, out float hu_ghost)
{
    float u_real = CalcU_Boundary(h_real, hu_real);
    float c_real = sqrt(boundary_g * h_real);
    float c_desired = sqrt(boundary_g * h_desired);
    float c_ghost = -u_real/2 - c_real + 2 * c_desired;

    if (c_ghost < 0) {
        h_ghost = 0;
        hu_ghost = hu_real + h_real * (2 * c_real - 4 * c_desired);
    } else {
        const float LIMIT = 2.0f;
        h_ghost = min(h_real + LIMIT, c_ghost*c_ghost / boundary_g);
        hu_ghost = 0;
    }
}    

float CalcSeaLevel( float2 tex_idx)
{
    return sea_level + (sa1 * cos(skx1 * tex_idx.x + sky1 * tex_idx.y - so1 * total_time)
                        + sa2 * cos(skx2 * tex_idx.x + sky2 * tex_idx.y - so2 * total_time)
                        + sa3 * cos(skx3 * tex_idx.x + sky3 * tex_idx.y - so3 * total_time)
                        + sa4 * cos(skx4 * tex_idx.x + sky4 * tex_idx.y - so4 * total_time)
                        ) * exp(-sdecay * tex_idx.y);
}

float4 NorthBoundary( VS_OUTPUT input ) : SV_TARGET
{
    int3 idx = int3(int(input.tex_idx.x), reflect_y - int(input.tex_idx.y), 0);
    float4 real = txState.Load(idx);
    float B = txBottom.Load(idx).b;

    float w_real = real.r;
    float hu_real = real.g;
    float hv_real = real.b;
    float h_real = w_real - B;

    float w_ghost;
    float hu_ghost;
    float hv_ghost;

    float SL = sea_level;
    
    if (idx.x >= inflow_x_min && idx.x <= inflow_x_max) {
        w_ghost = B + inflow_height;
        hu_ghost = 0;
        hv_ghost = inflow_height * (-inflow_speed);
        
    } else if (B > SL && solid_wall_flag) {
        w_ghost = w_real;
        hu_ghost = hu_real;
        hv_ghost = -hv_real;
    } else {
        FixedHBoundary(max(0, SL - B), h_real, hv_real, w_ghost, hv_ghost);
        w_ghost += B;
        hu_ghost = 0;
    }

    return float4(w_ghost, hu_ghost, hv_ghost, 0);
}
    
float4 EastBoundary( VS_OUTPUT input ) : SV_TARGET
{
    int3 idx = int3(reflect_x - int(input.tex_idx.x), int(input.tex_idx.y), 0);
    float4 real = txState.Load(idx);
    float B = txBottom.Load(idx).b;

    float w_real = real.r;
    float hu_real = real.g;
    float hv_real = real.b;
    float h_real = w_real - B;

    float w_ghost;
    float hu_ghost;
    float hv_ghost;

    float SL = CalcSeaLevel(input.tex_idx);
    
    if (B > SL && solid_wall_flag) {
        w_ghost = w_real;
        hu_ghost = -hu_real;
        hv_ghost = hv_real;
    } else {
        FixedHBoundary(max(0, SL - B), h_real, hu_real, w_ghost, hu_ghost);
        w_ghost += B;
        hv_ghost = 0;
    }

    return float4(w_ghost, hu_ghost, hv_ghost, 0);
}

float4 SouthBoundary( VS_OUTPUT input ) : SV_TARGET
{
    int3 idx = int3(int(input.tex_idx.x), 3 - int(input.tex_idx.y), 0);
    float4 real = txState.Load(idx);
    float B = txBottom.Load(idx).b;

    float w_real = real.r;
    float hu_real = real.g;
    float hv_real = real.b;
    float h_real = w_real - B;

    float w_ghost;
    float hu_ghost;
    float hv_ghost;

    float SL = CalcSeaLevel(input.tex_idx);    
    
    if (B > SL && solid_wall_flag) {
        w_ghost = w_real;
        hu_ghost = hu_real;
        hv_ghost = -hv_real;
    } else {
        FixedHBoundary(max(0, SL - B), h_real, -hv_real, w_ghost, hv_ghost);
        w_ghost += B;
        hv_ghost = -hv_ghost;
        hu_ghost = 0;
    }

    return float4(w_ghost, hu_ghost, hv_ghost, 0);
}    

float4 WestBoundary( VS_OUTPUT input ) : SV_TARGET
{
    // read inputs
    int3 idx = int3(3 - int(input.tex_idx.x), int(input.tex_idx.y), 0);
    float4 real = txState.Load(idx);
    float B = txBottom.Load(idx).b;
    
    float w_real = real.r;
    float hu_real = real.g;
    float hv_real = real.b;
    float h_real = w_real - B;
    
    float w_ghost;
    float hu_ghost;
    float hv_ghost;

    float SL = CalcSeaLevel(input.tex_idx);
    
    if (B > SL && solid_wall_flag) {
        w_ghost = w_real;
        hu_ghost = -hu_real;
        hv_ghost = hv_real;
    } else {
        FixedHBoundary(max(0, SL - B), h_real, -hu_real, w_ghost, hu_ghost);
        w_ghost += B;
        hu_ghost = -hu_ghost;
        hv_ghost = 0;
    }
    
    return float4(w_ghost, hu_ghost, hv_ghost, 0);
}






// GetStats shader
// This runs on a block of 4*4 pixels and returns aggregate stats
// back to the CPU
// The CPU then does the final aggregation over all the blocks.

// input: txState (t0), txBottom (t1)
// output: two RGBA float targets, containing the aggregate information.

struct GET_STATS_OUTPUT {
    float4 target0 : SV_TARGET0;    // sum(h),  sum(Bh + 0.5*h^2),  sum(h*u),  sum(h*v)
    float4 target1 : SV_TARGET1;    // sum(h*(u^2 + v^2)),  max(u^2+v^2),  max(h),  max((|u|+c)/dx, (|v|+c)/dy)
    float target2 : SV_TARGET2;    // max((u^2+v^2)/h)
};

GET_STATS_OUTPUT GetStats( VS_OUTPUT input )
{
    const int3 idx = int3(input.tex_idx.x, input.tex_idx.y, 0) * 4;

    float sum_h = 0;
    float sum_Bhh2 = 0;
    float sum_hu = 0;
    float sum_hv = 0;
    float sum_hu2v2 = 0;
    float max_u2v2 = 0;
    float max_h = 0;
    float max_cfl = 0;
    float max_f2 = 0;
    
    for (int j = 2; j < 6; ++j) {   // add 2 to avoid ghost zones
        for (int i = 2; i < 6; ++i) {
            const int3 idx2 = idx + int3(i, j, 0);
            
            const float4 state = txState.Load(idx2);
            const float w = state.r;
            const float hu = state.g;
            const float hv = state.b;

            const float3 B_here = txBottom.Load(idx2);
            const float B = B_here.b;

            const float h = max(0, w - B);
            const float c = sqrt(g * h);
            
            float u, v;
            float divide_by_h = CalcUV_Scalar(h, hu, hv, u, v);

            sum_h += h;
            sum_Bhh2 += h*(B + 0.5f * h);
            sum_hu += hu;
            sum_hv += hv;
            sum_hu2v2 += (hu * u + hv * v);

            float u2v2 = u*u + v*v;
            max_u2v2 = max(max_u2v2, u2v2);
            max_h = max(max_h, h);
            max_cfl = max(max_cfl, (abs(u) + c) * one_over_dx);
            max_cfl = max(max_cfl, (abs(v) + c) * one_over_dy);
            max_f2 = max(max_f2, u2v2 * divide_by_h);
        }
    }

    GET_STATS_OUTPUT output;
    output.target0 = float4(sum_h, sum_Bhh2, sum_hu, sum_hv);
    output.target1 = float4(sum_hu2v2, max_u2v2, max_h, max_cfl);
    output.target2 = max_f2;

    return output;
}
