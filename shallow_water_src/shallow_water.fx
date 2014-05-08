/* -*-c++-*-
 *
 * FILE:
 *   shallow_water.fx
 *
 * PURPOSE:
 *   Graphical shaders for shallow water demo
 *   (land, water, and sky)
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

// Textures

float4x4 gWorldViewProj : WorldViewProjection < string UIWidget = "None"; >; // world
float4x4 gWXf: World < string UIWidget = "None"; >; // world

// heightfield texture: 3 x 32-bit float
// r = B (height)
// g = dB/dx
// b = dB/dy

Texture2D<float3> txHeightfield : register( t0 );

// grass image texture.
Texture2D<float4> txGrass : register( t0 );

// a linear sampler
SamplerState samLinear : register( s0 );

// skybox cube texture
TextureCube txSkybox : register( t1 );


// water texture: 4 x 32-bit float
// r = w
// g = hu
// b = hv
// a = (unused)

Texture2D<float4> txWater : register( t1 );

// normal texture

// r = nX
// g = nY
// b = nZ
// a = unused

Texture2D<float4> txNormal : register( t2 );



// Constant Buffers

cbuffer MyConstBuffer : register( b0 )
{
    // used by renderer
    row_major float4x4 tex_to_clip;    // transforms (x_tex_idx, y_tex_idx, z_world, 1) to clip space
    float3 light_dir;   // in world space
    float ambient;
    float3 eye_mult, eye_trans;

    float2 terrain_tex_scale;

    // skybox matrix
    row_major float3x3 skybox_mtx;  // transforms world to clip space
    float4 pack;

    // refraction stuff
    float2 world_mult, world_trans;
    float2 world_to_grass_tex;
    float2 grass_tex_of_origin;

    // more settings
    float fresnel_coeff, fresnel_exponent;
    float specular_intensity, specular_exponent;
    float refractive_index;
    float attenuation_1, attenuation_2;
    
    int nx_plus_1, ny_plus_1;
    
    float3 deep_col;

    
};


// Structure definitions

struct VS_INPUT {
    int2 pos : INXX_POSITION;       // texture indices
};

struct TERRAIN_PS_INPUT {
    float4 pos : SV_POSITION;      // clip space
    float3 normal : NORMAL;        // world space
    float2 tex_coord : TEXCOORD;
};

struct WATER_PS_INPUT {
    float4 pos : SV_POSITION;      // clip space
    float3 normal : NORMAL;        // world space
    float3 eye : CAMERA_DIRECTION; // world space
    float water_depth : WATER_DEPTH;   // water depth (h) as float
    float2 world_pos : WORLD_XY_POS;   // world space
    float3 terrain_normal : TERRAIN_NORMAL;
};

struct SKYBOX_PS_INPUT {
    float4 pos : SV_POSITION;
    float3 tex : TEXCOORD;
};


// lighting function for terrain.

float3 TerrainColour(float3 tex_colour, float3 terrain_normal)
{
    float light = saturate(dot(light_dir, terrain_normal)) + ambient;
    return light * tex_colour;
}


// Vertex Shader for Terrain

TERRAIN_PS_INPUT TerrainVertexShader( VS_INPUT input )
{
    TERRAIN_PS_INPUT output;

    // lookup texture values at the input point
    const int3 tpos = int3(input.pos.x, input.pos.y, 0);
    const float3 tex_value = txHeightfield.Load( tpos );

    // compute position in clip space
    const float4 pos_in = float4( input.pos.x, input.pos.y, tex_value.r, 1 );
    output.pos = mul( tex_to_clip, pos_in );
    //output.pos = float4(output.pos.x * 0.5, output.pos.y * 0.1, output.pos.z, output.pos.a);

    const float3 norm = float3(-tex_value.g, -tex_value.b, 1.0f);
    output.normal = normalize(norm);

    // texture coords.
    // TODO: the tex coords might be better input in the vertex buffer rather than calculated here.
    output.tex_coord = terrain_tex_scale * float2(input.pos.x - 2, input.pos.y - 2);
    
    return output;
}


// Pixel Shader for Terrain

float4 TerrainPixelShader( TERRAIN_PS_INPUT input ) : SV_Target
{
    // tex lookup
    float3 tex_colour = txGrass.Sample( samLinear, input.tex_coord ).rgb;
    
    if(input.pos.x>100.0)
    {
        // Lighting calculation
        return float4(TerrainColour(tex_colour, input.normal), 1);
    }
    else
    {
        return float4(0.7, 0.2, 0.4, 1);
    }

}



// Vertex Shader for Water

WATER_PS_INPUT WaterVertexShader( VS_INPUT input )
{
    WATER_PS_INPUT output;

    // lookup the terrain level
    const int3 idx = int3(input.pos.x, input.pos.y, 0);
    const float3 ground_tex = txHeightfield.Load(idx);
    const float B = ground_tex.r;
    const float3 ground_normal = normalize(float3(-ground_tex.g, -ground_tex.b, 1));
    
    // lookup the water level
    float w = txWater.Load(idx).r;

    // simple way to ensure there is no "gap" at the edge of the box
    if (idx.x==2 || idx.y==2 || idx.x==nx_plus_1 || idx.y==ny_plus_1) w=B;

    // calculate water depth, this is used in the pixel shader
    output.water_depth = (w - B);

    // compute position in clip space
    const float dmin = 0.05f;
    const float vert_bias = (output.water_depth < dmin ? 2*(output.water_depth - dmin) : 0);
    const float4 pos_in = float4( input.pos.x, input.pos.y, w + vert_bias, 1 );
    output.pos = mul( tex_to_clip, pos_in );
    //const float4 pos_in = float4( input.pos.x*0.3, input.pos.y*0.3, w + vert_bias, 1 );
    //output.pos = mul(gWXf, pos_in);


    // now compute the normal
    output.normal = txNormal.Load(idx).rgb;

    // now compute the eye vector (from vertex towards camera)
    output.eye = normalize(eye_mult * pos_in.xyz + eye_trans);
    
    // send the world pos (xy) through, for refraction calculations
    output.world_pos = world_mult * pos_in.xy + world_trans;

    output.terrain_normal = ground_normal;
    
    return output;
}



// simplified gamma correction.

float3 ToLinear(float3 srgb)
{
    return pow(abs(srgb), 2.2f);
}

float3 FromLinear(float3 lin)
{
    return pow(abs(lin), 1.0f / 2.2f);
}




// Pixel Shader for Water

float4 WaterPixelShader( WATER_PS_INPUT input ) : SV_Target
{
    // approximate Fresnel factor
    // (this is the percentage that reflects, as opposed to transmits)
    float fresnel = (1-fresnel_coeff) + fresnel_coeff * pow(max(0, 1 - dot(input.eye, input.normal)), fresnel_exponent);
    
    // reflected light
    float3 reflection_dir = reflect(-input.eye, input.normal);
    float3 reflect_col = txSkybox.Sample(samLinear, reflection_dir).rgb;

    //specular
    reflect_col += specular_intensity * pow(max(0,dot(reflection_dir, light_dir)), specular_exponent);

    // refracted light
    float3 refract_col;
    float attenuation_factor;
    float3 refract_dir = refract( -input.eye, input.normal, 1.0f / refractive_index);

    if (refract_dir.z > 0) {
        // refracted ray goes up into the sky
        // (this can only happen with slanted water surfaces)
        // note: we have no way of knowing the distance through the water that the
        // ray has travelled, so we just hard code the attenuation
        refract_col = txSkybox.Sample(samLinear, refract_dir).rgb;
        attenuation_factor = 0.8f;
        
    } else {
        // refracted ray hits the grass
        float dist = -input.water_depth / refract_dir.z;
        float2 base_pos = input.world_pos.xy + dist * refract_dir.xy;
        float2 refract_tex_coord = world_to_grass_tex * base_pos + grass_tex_of_origin;
        refract_col = txGrass.Sample(samLinear, refract_tex_coord).rgb;
        refract_col = TerrainColour(refract_col, input.terrain_normal);  // apply terrain lighting (approximation)
        attenuation_factor = (1-attenuation_1) * exp(-attenuation_2 * dist);
    }
    
    // combine reflection & refraction
    float3 result = FromLinear(lerp(lerp(ToLinear(deep_col),
                                         ToLinear(refract_col),
                                         attenuation_factor),
                                    ToLinear(reflect_col),
                                    fresnel));
    
    return float4(result, 1);
}


// Skybox shaders

SKYBOX_PS_INPUT SkyboxVertexShader( float3 pos : POSITION )
{
    SKYBOX_PS_INPUT output;

    output.pos = mul(skybox_mtx, pos).xyzz;   // set w=z
    output.pos.w *= 1.0001f;  // prevent it colliding with the far clip plane

    output.tex = pos;

    return output;
}

float4 SkyboxPixelShader( SKYBOX_PS_INPUT input ) : SV_TARGET
{
    // tex lookup
    return txSkybox.Sample(samLinear, input.tex);
}

struct VS_OUTPUT {
    float4 pos : SV_POSITION;
    float2 tex_idx : TEX_IDX;  // should be cast to an integer by the pixel shader
    float radius : RADIUS;
};

VS_OUTPUT main(VS_INPUT vs_in)
 {
   VS_OUTPUT vs_out = (VS_OUTPUT)0;
   
   return vs_out;
 } 