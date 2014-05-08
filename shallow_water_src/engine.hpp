/*
 * FILE:
 *   engine.hpp
 *
 * PURPOSE:
 *   "Engine" class for the shallow water demo. Manages all the
 *   Direct3D objects and contains top-level rendering routines.
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

#ifndef ENGINE_HPP
#define ENGINE_HPP

#include "settings.hpp"

#include "coercri/dx11/core/com_ptr_wrapper.hpp"

#include <d3d11.h>
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

class ShallowWaterEngine {
public:
    ShallowWaterEngine(ID3D11Device *device_, ID3D11DeviceContext *context_);

    void newTerrainSettings();  // call if any terrain settings change.
    void remesh(ResetType rt);   // call if mesh size changes.

    // call following camera move or window resize
    void moveCamera(float cam_x, float cam_y, float cam_z, float yaw, float pitch, int vp_width, int vp_height);

    // update the water
    void timestep();

    // set timestep to the given dt (multiplied by time_acceleration),
    // or to safety_factor * CFL-timestep,
    // whichever is smaller.
    void resetTimestep(float dt);
    
    // render a frame to the given render target
    void render(ID3D11RenderTargetView *rtv);


    // mouse picking
    bool mousePick(int mouse_x, int mouse_y,
                   float &world_x, float &world_y, float &world_z, float &depth, float &u, float &v);
    void applyMouseShader(float world_x, float world_y, float dt);

    // get water height (h + B) at the given point
    float getWaterHeight(float world_x, float world_y);
    
    
private:
    void createShadersAndInputLayout();
    void createMeshBuffers();
    void createSimBuffers();
    
    void createTerrainTexture();
    void fillTerrainTexture();
    void fillTerrainTextureLite();
    void createSimTextures(ResetType reset_type);
    void createConstantBuffers();
    void fillConstantBuffers();
    void createDepthStencil(int w, int h);
    void loadGraphics();
    void loadSkybox();

    void setupMousePicking();
    void raiseLowerTerrain(float wx, float wy, float dt);
    
private:
    ID3D11Device *device;
    ID3D11DeviceContext *context;

    // vertex & pixel shaders
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psTerrainVertexShader, m_psWaterVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psTerrainPixelShader, m_psWaterPixelShader;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psSimVertexShader;  // shared between KP07 and Lax-Wendroff methods
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psSimPixelShader[3], m_psGetStatsPixelShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psLaxWendroffPixelShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psBoundaryPixelShader[4];

    // input layouts
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psInputLayout, m_psSimInputLayout;
    
    // vertex and index buffers (for the triangle meshes).
    // can be used for both terrain & water.
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psMeshVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psMeshIndexBuffer;

    // vertex buffers for the simulation render-to-texture.
    // (one for each pass.)
    // TODO: Might be better to have one large vertex buffer, with offsets,
    // rather than lots of little ones like this.
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psSimVertexBuffer11, m_psSimVertexBuffer10, m_psSimVertexBuffer00;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psGetStatsVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psGetStatsIndexBuffer;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psBoundaryVertexBuffer[4];

    // constant buffers
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psConstantBuffer, m_psSimConstantBuffer, m_psBoundaryConstantBuffer;
    
    // textures:
    //  -- terrain (contains B, dB/dx, dB/dy)
    //  -- bottom (contains BY, BX, BA; used for simulation)
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psTerrainTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psTerrainTextureView;
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psBottomTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psBottomTextureView;

    // simulation textures:
    // [sim_idx] = state
    // [1-sim_idx] = output state / H
    // [2] = U
    // [3] = V
    // [4] = XFLUX
    // [5] = YFLUX
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psSimTexture[7], m_psGetStatsTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psSimTextureView[7];
    Coercri::ComPtrWrapper<ID3D11RenderTargetView> m_psSimRenderTargetView[7], m_psGetStatsRenderTargetView;
    int sim_idx;
    bool bootstrap_needed;

    // staging textures
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psFullSizeStagingTexture;
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psGetStatsStagingTexture4, m_psGetStatsStagingTexture1;

    // graphical textures
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psGrassTexture;
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psGrassTextureView;
    Coercri::ComPtrWrapper<ID3D11SamplerState> m_psLinearSamplerState;

    // skybox
    Coercri::ComPtrWrapper<ID3D11ShaderResourceView> m_psSkyboxView;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psSkyboxVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psSkyboxPixelShader;
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psSkyboxInputLayout;
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psSkyboxVertexBuffer;
    
    // depth buffer stuff
    Coercri::ComPtrWrapper<ID3D11Texture2D> m_psDepthStencil;
    Coercri::ComPtrWrapper<ID3D11DepthStencilView> m_psDepthStencilView;

    // mouse picking
    Coercri::ComPtrWrapper<ID3D11Buffer> m_psLeftMouseConstantBuffer, m_psLeftMouseVertexBuffer;
    Coercri::ComPtrWrapper<ID3D11VertexShader> m_psLeftMouseVertexShader;
    Coercri::ComPtrWrapper<ID3D11PixelShader> m_psLeftMousePixelShader;
    Coercri::ComPtrWrapper<ID3D11InputLayout> m_psLeftMouseInputLayout;
    
    // current timestep
    float current_timestep;
    float total_time;

    // viewport/camera state
    int vp_width, vp_height;
    float pitch, yaw, camera_x, camera_y, camera_z;
    float xpersp, ypersp;  // perspective coefficients.
    float near_plane_dist, far_plane_dist;
};

#endif
