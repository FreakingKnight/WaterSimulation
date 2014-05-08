/*
 * FILE:
 *   engine.cpp
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
 * NOTES:
 *   Most of the code for the demo ended up in this class
 *   (ShallowWaterEngine). That's probably okay in this case, since
 *   this is only a "prototype", but for a real implementation we
 *   might consider refactoring / splitting up into different classes
 *   perhaps.
 *   
 */

#include "engine.hpp"
#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include "coercri/dx11/core/dx_error.hpp"
#include "coercri/gfx/load_bmp.hpp"
#include "coercri/gfx/pixel_array.hpp"

#include "boost/scoped_array.hpp"

#include <d3d11.h>
#include <d3dx11.h>
#include <d3dcompiler.h>
#include <d3dx10math.h>
#include <xnamath.h>

#include <cmath>
#include <fstream>
#include <string>

#include <d3d9.h>

#ifdef near
#undef near
#endif
#ifdef far
#undef far
#endif

// Turning off USE_KP07 activates an experimental Lax-Wendroff solver.
// Unfortunately this is buggy and unstable currently, so leaving USE_KP07 
// on is recommended.
#define USE_KP07

// Debugging switches
//#define DUMP_TO_FILE
//#define RUN_CHECKS

#if defined(DUMP_TO_FILE) || defined(RUN_CHECKS)
#include <iomanip>
#endif

namespace {

    const float PI = 4.0f * std::atan(1.0f);

    float GaussianBrush(float x, float y, float r)
    {
        const float sigma = 1.0f / 3.0f;
        float z = (x*x + y*y) / (r*r);
        if (z > 1) return 0;
        else return 2 * std::exp(-z/(2*sigma*sigma));
    }
    
    class MapTexture {
    public:
        MapTexture(ID3D11DeviceContext &cxt, ID3D11Texture2D & tex)
            : context(cxt), texture(tex)
        {
            HRESULT hr = context.Map(&texture, 0, D3D11_MAP_READ, 0, &msr);
            if (FAILED(hr)) {
                throw Coercri::DXError("Map failed", hr);
            }
        }

        ~MapTexture()
        {
            context.Unmap(&texture, 0);
        }

        // allow access to the pointer / pitch values
        D3D11_MAPPED_SUBRESOURCE msr;

    private:
        ID3D11DeviceContext &context;
        ID3D11Texture2D &texture;  // must be a staging texture
    };                   
    
    // VS input for water/terrain rendering
    struct MeshVertex {
        int x, y;
    };

    // Const buffer for water/terrain rendering
    struct MyConstBuffer {
        XMMATRIX tex_to_clip;   // transforms (tex_x, tex_y, world_z, 1) into clip space 
        XMFLOAT3 light_dir;         // in world space 
        float ambient;
        XMFLOAT3 eye_mult;
        float pack6;
        XMFLOAT3 eye_trans;
        float pack7;
        XMFLOAT2 terrain_tex_scale;
        float pack8, pack9;
        XMMATRIX skybox_mtx;
        float world_mult_x, world_mult_y, world_trans_x, world_trans_y;
        float world_to_grass_tex_x, world_to_grass_tex_y;
        float grass_tex_of_origin_x, grass_tex_of_origin_y;

        float fresnel_coeff, fresnel_exponent;
        float specular_intensity, specular_exponent;
        float refractive_index;
        float attenuation_1, attenuation_2;
        
        int nx_plus_1, ny_plus_1;

        float deep_r, deep_g, deep_b;
        
    };    

    // Const buffer used for the simulation
    struct SimConstBuffer {
        float two_theta;
        float two_over_nx_plus_four;
        float two_over_ny_plus_four;
        float g;
        float half_g;
        float g_over_dx;
        float g_over_dy;
        float one_over_dx;
        float one_over_dy;
        float dt;
        float epsilon;      // usually dx^4
        int nx, ny;
        float friction;  // m s^-1
    };

    // Const buffer used for boundary conditions
    struct BoundaryConstBuffer {
        float boundary_epsilon;
        int reflect_x, reflect_y;
        int solid_wall_flag;
        int inflow_x_min, inflow_x_max;
        float sea_level, inflow_height, inflow_speed;
        float g;
        float total_time;
        float sa1, skx1, sky1, so1;
        float sa2, skx2, sky2, so2;
        float sa3, skx3, sky3, so3;
        float sa4, skx4, sky4, so4;
        float sdecay;
    };

    struct LeftMouseConstBuffer {
        float scale_x, scale_y;
        float bias_x, bias_y;
        float two_over_nx_plus_four;
        float two_over_ny_plus_four;
        float disp_A, disp_B;
    };

    float CalcEpsilon()
    {
        const float W = GetSetting("valley_width");
        const float L = GetSetting("valley_length");
        const float nx = GetSetting("mesh_size_x");
        const float ny = GetSetting("mesh_size_y");
        const float dx = W / (nx-1);
        const float dy = L / (ny-1);
        return std::min(0.5f, std::pow(std::max(dx, dy), 4));
    }

    float CalcU(float h, float hu)
    {
        float epsilon = CalcEpsilon();
        float h2 = h*h;
        float h4 = h2*h2;
        float divide_by_h = sqrt(2.0f) * h / sqrt(h4 + std::max(h4, epsilon));
        return divide_by_h * hu;
    }
    
    void CompileShader(const char *filename,
                       const char *entry_point,
                       const char *profile,
                       Coercri::ComPtrWrapper<ID3DBlob> &compiled_code)
    {
        DWORD compile_flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_WARNINGS_ARE_ERRORS;
#ifdef _DEBUG
        compile_flags |= D3DCOMPILE_DEBUG;
#endif

        ID3DBlob *output_blob = 0;
        ID3DBlob *error_blob = 0;
        
        HRESULT hr = D3DX11CompileFromFile(filename,
                                           0,   // defines
                                           0,   // include handler
                                           entry_point,
                                           profile,
                                           compile_flags,
                                           0,  // effect flags (ignored)
                                           0,  // thread pump
                                           &output_blob,
                                           &error_blob,
                                           0);    // used for async compilation only

        compiled_code.reset(output_blob);
        Coercri::ComPtrWrapper<ID3DBlob> error_sentinel(error_blob);
        
        if (FAILED(hr)) {
            std::string msg = "D3DX11CompileFromFile failed";
            if (error_blob) {
                OutputDebugStringA((char*)error_blob->GetBufferPointer());
                msg += "\n\n";
                msg += (char*)error_blob->GetBufferPointer();
            }
            throw Coercri::DXError(msg, hr);
        }
    }

    void CreateVertexShader(ID3D11Device *device,
                            const char *filename,
                            const char *entry_point,
                            Coercri::ComPtrWrapper<ID3DBlob> &bytecode,
                            Coercri::ComPtrWrapper<ID3D11VertexShader> &output)
    {
        const char * vs_profile = "vs_4_0";

        CompileShader(filename, entry_point, vs_profile, bytecode);

        ID3D11VertexShader *vert_shader = 0;
        HRESULT hr = device->CreateVertexShader(bytecode->GetBufferPointer(),
                                                bytecode->GetBufferSize(),
                                                0,
                                                &vert_shader);
        if (FAILED(hr)) {
            throw Coercri::DXError("CreateVertexShader failed", hr);
        }
        output.reset(vert_shader);
    }

    void CreatePixelShader(ID3D11Device *device,
                           const char *filename,
                           const char *entry_point,
                           Coercri::ComPtrWrapper<ID3D11PixelShader> &output)
    {
        const char * ps_profile = "ps_4_0";

        Coercri::ComPtrWrapper<ID3DBlob> bytecode;
        CompileShader(filename, entry_point, ps_profile, bytecode);

        ID3D11PixelShader *pixel_shader = 0;
        HRESULT hr = device->CreatePixelShader(bytecode->GetBufferPointer(),
                                               bytecode->GetBufferSize(),
                                               0,
                                               &pixel_shader);
        if (FAILED(hr)) {
            throw Coercri::DXError("CreatePixelShader failed", hr);
        }
        output.reset(pixel_shader);
    }

    void GetTextureSize(ID3D11Texture2D *tex, int &width, int &height)
    {
        D3D11_TEXTURE2D_DESC td;
        tex->GetDesc(&td);

        width = td.Width;
        height = td.Height;
    }

    void GetRenderTargetSize(ID3D11RenderTargetView *view, int &width, int &height)
    {
        ID3D11Resource *resource;
        view->GetResource(&resource);
        Coercri::ComPtrWrapper<ID3D11Resource> psResource(resource);

        ID3D11Texture2D *texture;
        resource->QueryInterface(__uuidof(ID3D11Texture2D), reinterpret_cast<void**>(&texture));
        Coercri::ComPtrWrapper<ID3D11Texture2D> psTexture(texture);

        GetTextureSize(texture, width, height);
    }

    void CreateTextureImpl(ID3D11Device *device,
                           const D3D11_TEXTURE2D_DESC &td,
                           const D3D11_SUBRESOURCE_DATA *srd,
                           Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                           Coercri::ComPtrWrapper<ID3D11ShaderResourceView> *out_srv,
                           Coercri::ComPtrWrapper<ID3D11RenderTargetView> *out_rtv)
    {
        ID3D11Texture2D *pTexture;
        HRESULT hr = device->CreateTexture2D(&td, srd, &pTexture);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create texture", hr);
        }
        out_tex.reset(pTexture);

        if (out_srv) {
            D3D11_SHADER_RESOURCE_VIEW_DESC sd;
            memset(&sd, 0, sizeof(sd));
            sd.Format = td.Format;
            sd.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
            sd.Texture2D.MipLevels = td.MipLevels;
            
            ID3D11ShaderResourceView *pSRV;
            hr = device->CreateShaderResourceView(pTexture, &sd, &pSRV);
            if (FAILED(hr)) {
                throw Coercri::DXError("Failed to create shader resource view for texture", hr);
            }
            out_srv->reset(pSRV);
        }

        if (out_rtv) {
            D3D11_RENDER_TARGET_VIEW_DESC rd;
            memset(&rd, 0, sizeof(rd));
            rd.Format = td.Format;
            rd.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
            
            ID3D11RenderTargetView *rtv;
            HRESULT hr = device->CreateRenderTargetView(pTexture, &rd, &rtv);
            if (FAILED(hr)) {
                throw Coercri::DXError("Failed to create render target view", hr);
            }
            out_rtv->reset(rtv);
        }
    }
    

    // Creates an texture, with optional initial data,
    // and creates optional shader resource view, and optional render target view
    void CreateTexture(ID3D11Device * device,
                       int width,
                       int height,
                       const D3D11_SUBRESOURCE_DATA *initial_data,
                       DXGI_FORMAT format,
                       bool staging,
                       Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                       Coercri::ComPtrWrapper<ID3D11ShaderResourceView> *out_srv,
                       Coercri::ComPtrWrapper<ID3D11RenderTargetView> *out_rtv)
    {
        D3D11_TEXTURE2D_DESC td;
        memset(&td, 0, sizeof(td));
        td.Width = width;
        td.Height = height;
        td.MipLevels = 1;
        td.ArraySize = 1;
        td.Format = format;
        td.SampleDesc.Count = 1;
        td.Usage = D3D11_USAGE_DEFAULT;

        if (out_srv) td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        if (out_rtv) td.BindFlags |= D3D11_BIND_RENDER_TARGET;

        if (staging) {
            td.Usage = D3D11_USAGE_STAGING;
            td.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        }

        CreateTextureImpl(device, td, initial_data, out_tex, out_srv, out_rtv);
    }

    int GetNumMipLevels(int width, int height)
    {
        int num_levels = 1;
        while (width > 1 && height > 1) {
            width /= 2;
            height /= 2;
            ++num_levels;
        }
        return num_levels;
    }

    int TexIndex(int x, int y, int component, int width)
    {
        return (y * width + x)*4 + component;
    }
    
    // creates an IMMUTABLE texture and auto generates the mip maps.
    void CreateTextureWithMips(ID3D11Device *device,
                               int width,
                               int height,
                               const unsigned char * initial_data,  // 32-bit RGBA (8-bit per channel)
                               Coercri::ComPtrWrapper<ID3D11Texture2D> &out_tex,
                               Coercri::ComPtrWrapper<ID3D11ShaderResourceView> &out_view)
    {
        const int num_mip_levels = GetNumMipLevels(width, height);

        D3D11_TEXTURE2D_DESC td;
        memset(&td, 0, sizeof(td));
        td.Width = width;
        td.Height = height;
        td.MipLevels = num_mip_levels;
        td.ArraySize = 1;
        td.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        td.SampleDesc.Count = 1;
        td.Usage = D3D11_USAGE_IMMUTABLE;
        td.BindFlags = D3D11_BIND_SHADER_RESOURCE;

        boost::scoped_array<D3D11_SUBRESOURCE_DATA> srd(new D3D11_SUBRESOURCE_DATA[num_mip_levels]);
        boost::scoped_array<std::vector<unsigned char> > pixels(new std::vector<unsigned char>[num_mip_levels]);

        srd[0].pSysMem = initial_data;
        srd[0].SysMemPitch = width * 4;

        for (int level = 1; level < num_mip_levels; ++level) {

            const int prev_width = width;
            width = std::max(width/2, 1);
            height = std::max(height/2, 1);

            pixels[level].resize(width*height*4);
            unsigned char * out = &pixels[level][0];
            
            const unsigned char * in;
            if (level == 1) {
                in = initial_data;
            } else {
                in = &pixels[level-1][0];
            }
            
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < height; ++x) {
                    for (int component = 0; component < 4; ++component) {

                        const int in1 = in[TexIndex(2*x,   2*y,   component, prev_width)];
                        const int in2 = in[TexIndex(2*x+1, 2*y,   component, prev_width)];
                        const int in3 = in[TexIndex(2*x,   2*y+1, component, prev_width)];
                        const int in4 = in[TexIndex(2*x+1, 2*y+1, component, prev_width)];

                        const int avg = (in1 + in2 + in3 + in4 + 2) / 4;

                        *out++ = avg;
                    }
                }
            }

            srd[level].pSysMem = &pixels[level][0];
            srd[level].SysMemPitch = width * 4;
        }

        CreateTextureImpl(device, td, &srd[0], out_tex, &out_view, 0);
    }
    
    void CreateSimBuffer(ID3D11Device *device, Coercri::ComPtrWrapper<ID3D11Buffer> &vert_buf, 
                         int i_left, int i_top, int i_right, int i_bottom)
    {
        // assumes viewport of (nx+4) * (ny+4),
        // and renders all pixels from (left,top) (inclusive) to (right,bottom) (exclusive).
         
        // the tex coords sent to the pixel shader are e.g. (0.5f, 0.5f) for the top left pixel in the render target.
        
        const float left = float(i_left);
        const float right = float(i_right);
        const float top = float(i_top);
        const float bottom = float(i_bottom);

        const float vertices[12] = {
            left, top,
            right, top,
            left, bottom,
            right, top,
            right, bottom,
            left, bottom
        };
                
        // create the vertex buffer
        D3D11_BUFFER_DESC bd;
        memset(&bd, 0, sizeof(bd));
        bd.ByteWidth = 12 * sizeof(float);
        bd.Usage = D3D11_USAGE_IMMUTABLE;
        bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

        D3D11_SUBRESOURCE_DATA sd;
        memset(&sd, 0, sizeof(sd));
        sd.pSysMem = &vertices[0];
        
        ID3D11Buffer *pBuffer;
        HRESULT hr = device->CreateBuffer(&bd, &sd, &pBuffer);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create vertex buffer for the simulation", hr);
        }
        vert_buf.reset(pBuffer);
    }

    void CreateSimBuffer( ID3D11Device *device, Coercri::ComPtrWrapper<ID3D11Buffer> &vert_buf, 
                          Coercri::ComPtrWrapper<ID3D11Buffer> &indices_buf, 
                          int i_left, int i_top, int i_right, int i_bottom)
    {
        // assumes viewport of (nx+4) * (ny+4),
        // and renders all pixels from (left,top) (inclusive) to (right,bottom) (exclusive).
         
        // the tex coords sent to the pixel shader are e.g. (0.5f, 0.5f) for the top left pixel in the render target.
        
        const float left = float(i_left);
        const float right = float(i_right);
        const float top = float(i_top);
        const float bottom = float(i_bottom);

        const float vertices[12] = {
            left, top,
            right, top,
            left, bottom,
            right, top,
            right, bottom,
            left, bottom
        };
                
        // create the vertex buffer
        D3D11_BUFFER_DESC bd;
        memset(&bd, 0, sizeof(bd));
        bd.ByteWidth = 12 * sizeof(float);
        bd.Usage = D3D11_USAGE_IMMUTABLE;
        bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

        D3D11_SUBRESOURCE_DATA sd;
        memset(&sd, 0, sizeof(sd));
        sd.pSysMem = &vertices[0];
        
        ID3D11Buffer *pBuffer;
        HRESULT hr = device->CreateBuffer(&bd, &sd, &pBuffer);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create vertex buffer for the simulation", hr);
        }
        vert_buf.reset(pBuffer);

        // create data for the index buffer
        boost::scoped_array<int> indices(new int[6]);
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        indices[3] = 3;
        indices[4] = 4;
        indices[5] = 5;
            
        // create the index buffer
        bd.ByteWidth = sizeof(int) * 6;
        bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
        sd.pSysMem = &indices[0];
        sd.SysMemPitch = 0;

        hr = device->CreateBuffer(&bd, &sd, &pBuffer);
        if (FAILED(hr)) {
            throw Coercri::DXError("Failed to create mesh index buffer", hr);
        }
        indices_buf.reset(pBuffer);

    }


    int RoundUpTo16(int x)
    {
        return (x + 15) & (~15);
    }

#ifdef RUN_CHECKS

    std::ofstream g_chk("c:/cygwin/home/stephen/projects/shallow-water/runtime_checks.txt");

    void RunChecks(ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex)
    {
        static int count = 0;
        
        const int nx = GetIntSetting("mesh_size_x");
        const int ny = GetIntSetting("mesh_size_y");
        
        if (count == 0) {
            // check the BX, BY textures average to BA
            float max_err = 0;
            for (int j = 2; j < ny+2; ++j) {
                for (int i = 2; i < nx+2; ++i) {
                    const BottomEntry &b = g_bottom[j*(nx+4)+i];
                    const BottomEntry &bw = g_bottom[j*(nx+4)+i-1];
                    const BottomEntry &bs = g_bottom[(j-1)*(nx+4)+i];
                    
                    max_err = std::max(max_err, std::fabs(
                        b.BA - 0.25f*(b.BX + b.BY + bw.BX + bs.BY)));
                }
            }
            g_chk << "B err = " << max_err << "\n";
        }

        context->CopyResource(staging, tex);

        MapTexture m(staging);

        float h_min = 1e10;
        int h_min_i, h_min_j;
        float h_tot = 0;
        float w_tot = 0;
        float u_max = -1e10;  // max(|u|)
        float v_max = -1e10;  // max(|v|)

        bool bad_flag = false;
        
        for (int j = 2; j < ny + 2; ++j) {
            for (int i = 2; i < nx + 2; ++i) {
                const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
                const float *p = reinterpret_cast<const float*>(q) + i * 4;

                const float w = p[0];
                const float hu = p[1];
                const float hv = p[2];

                const float h = w - g_bottom[(nx+4)*j + i].BA;
                const float h2 = h*h;
                const float h4 = h2*h2;
                const float div_by_h = sqrt(2.0f) * h / sqrt(h4 + std::max(h4, 0.06f)); // epsilon = 0.06 hard coded !!
                const float u = hu * div_by_h;
                const float v = hv * div_by_h;

                if (h < h_min) {
                    h_min = h;
                    h_min_i = i;
                    h_min_j = j;
                }
                u_max = std::max(std::fabs(u), u_max);
                v_max = std::max(std::fabs(v), v_max);
                h_tot += h;
                w_tot += w;

                if (!_finite(w) || !_finite(hu) || !_finite(hv)) bad_flag = true;
            }
        }

        g_chk << std::setw(6)  << count++ 
              << std::setw(14) << h_min 
              << std::setw(5) << h_min_i
              << std::setw(5) << h_min_j
              << std::setw(14) << u_max 
              << std::setw(14) << v_max 
              << std::setw(14) << h_tot 
              << std::setw(14) << w_tot << "\n";
        if (bad_flag) {
            g_chk << "NON FINITE VALUE DETECTED, STOPPING\n";
            g_chk.close();
            std::abort();
        }
    }
#endif
    
#ifdef DUMP_TO_FILE
    void DumpToFile(std::ofstream &str, ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex)
    {
        const int nx = GetIntSetting("mesh_size_x");
        const int ny = GetIntSetting("mesh_size_y");
        
        context->CopyResource(staging, tex);

        MapTexture m(staging);
        
        for (int j = 0; j < ny + 4; ++j) {
            for (int i = 0; i < nx + 4; ++i) {
                const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
                const float *p = reinterpret_cast<const float*>(q) + i * 4;
                str << i-2 << "\t" << j-2 << "\t" << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << p[3] << "\n";
            }
        }
    }

    float GetSum(ID3D11DeviceContext *context, ID3D11Texture2D *staging,
                 int xmin, int xmax, int ymin, int ymax)
    {
        float result = 0;
        
        const int nx = GetIntSetting("mesh_size_x");
        const int ny = GetIntSetting("mesh_size_y");

        MapTexture m(staging);

        for (int j = ymin+2; j <= ymax+2; ++j) {
            for (int i = xmin+2; i <= xmax+2; ++i) {
                const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
                const float *p = reinterpret_cast<const float*>(q) + i * 4;
                result += *p;
            }
        }

        return result;
    }

    void WriteTotal(std::ofstream &str, ID3D11DeviceContext *context, ID3D11Texture2D *staging, ID3D11Texture2D *tex)
    {
        const int nx = GetIntSetting("mesh_size_x");
        const int ny = GetIntSetting("mesh_size_y");
        
        context->CopyResource(staging, tex);

        MapTexture m(staging);

        float total = 0;

        for (int j = 0; j < ny + 4; ++j) {
            for (int i = 0; i < nx + 4; ++i) {
                const char *q = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
                const float *p = reinterpret_cast<const float*>(q) + i * 4;
                total += *p;
            }
        }

        str << "Total H = " << total << "\n\n";
    }
#endif

    void SeaWaveSettings(char c, float &sa, float &skx, float &sky, float &so)
    {
        using std::string;
        sa = GetSetting(string("sa") + c);
        float k = GetSetting(string("sk") + c);
        float kdir = GetSetting(string("sk") + c + "_dir");
        skx = k * std::cos(kdir);
        sky = k * std::sin(kdir);
        so = GetSetting(string("so") + c);
    }
}

ShallowWaterEngine::ShallowWaterEngine(ID3D11Device *device_,
                                       ID3D11DeviceContext *context_)
    : device(device_), context(context_), current_timestep(0), total_time(0)
{
    // create D3D objects
    createShadersAndInputLayout();
    createConstantBuffers();
    loadGraphics();
    //createBlendState();
    loadSkybox();
    setupMousePicking();
    
    remesh(R_MESH);
    newTerrainSettings();
    moveCamera(0, 0, 30, 0, 0, 100, 100);
}

void ShallowWaterEngine::remesh(ResetType reset_type)
{
    createMeshBuffers();
    createSimBuffers();
    createTerrainTexture();
    fillTerrainTextureLite();
    createSimTextures(reset_type);
}

void ShallowWaterEngine::newTerrainSettings()
{
    fillTerrainTexture();
}

void ShallowWaterEngine::moveCamera(float x, float y, float z, float yaw_, float pitch_, int vw, int vh)
{
    pitch = pitch_;
    yaw = yaw_;
    camera_x = x;
    camera_y = y;
    camera_z = z;
    
    vp_width = vw;
    vp_height = vh;

    near_plane_dist = 1;
    far_plane_dist = 2 * std::max(GetSetting("valley_width"), GetSetting("valley_length"));

    const float fov = GetSetting("fov");  // x fov in degrees
    const float fov_r = fov / 45.0f * std::atan(1.0f);
    xpersp = 1.0f / tan(fov_r * 0.5f);
    ypersp = float(vw)/float(vh) * xpersp;
    //ypersp = 1.0f / tan(fov_r * 0.5f);
    //xpersp = float(vw)/float(vh) * ypersp;
    
    fillConstantBuffers();
}    

void ShallowWaterEngine::timestep()
{
    D3DPERF_BeginEvent(0, L"demo_timestep");

    /// Allow only a single timestep for debugging
    //static bool debug_flag = false;
    //if (debug_flag) return;
    //debug_flag = true;
    ///////////////////////////////////

    // Get some resource pointers
    
    ID3D11Buffer * cst_buf = m_psSimConstantBuffer.get();

    ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
    ID3D11ShaderResourceView * old_state_tex = m_psSimTextureView[sim_idx].get();
    ID3D11ShaderResourceView * new_state_or_h_tex = m_psSimTextureView[1 - sim_idx].get();
    ID3D11ShaderResourceView * u_tex = m_psSimTextureView[2].get();
    ID3D11ShaderResourceView * v_tex = m_psSimTextureView[3].get();
    ID3D11ShaderResourceView * xflux_tex = m_psSimTextureView[4].get();
    ID3D11ShaderResourceView * yflux_tex = m_psSimTextureView[5].get();

    ID3D11RenderTargetView * old_state_or_h_target = m_psSimRenderTargetView[sim_idx].get();
    ID3D11RenderTargetView * new_state_or_h_target = m_psSimRenderTargetView[1 - sim_idx].get();
    ID3D11RenderTargetView * u_target = m_psSimRenderTargetView[2].get();
    ID3D11RenderTargetView * v_target = m_psSimRenderTargetView[3].get();
    ID3D11RenderTargetView * xflux_target = m_psSimRenderTargetView[4].get();
    ID3D11RenderTargetView * yflux_target = m_psSimRenderTargetView[5].get();
    ID3D11RenderTargetView * normal_target = m_psSimRenderTargetView[6].get();

    
    // Common Settings
    
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psSimInputLayout.get());

    context->VSSetShader(m_psSimVertexShader.get(), 0, 0);
    
    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    vp.Width = float(nx + 4);
    vp.Height = float(ny + 4);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    context->VSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetConstantBuffers(0, 1, &cst_buf);

    ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

    ID3D11ShaderResourceView *pNULL = 0;

#ifdef USE_KP07
    if (bootstrap_needed) {
        // Pass 1
        // read: old_state; bottom
        // write: h, u, v
        
        vert_buf = m_psSimVertexBuffer11.get();
        context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
        
        ID3D11RenderTargetView * p1_tgt[] = { new_state_or_h_target, u_target, v_target, normal_target };
        context->OMSetRenderTargets(4, &p1_tgt[0], 0);
        
        context->PSSetShader(m_psSimPixelShader[0].get(), 0, 0);
        context->PSSetShaderResources(0, 1, &old_state_tex);
        context->PSSetShaderResources(1, 1, &bottom_tex);
        
        context->Draw(6, 0);

        bootstrap_needed = false;
    }

    
#ifdef DUMP_TO_FILE

    // delete the file first time
    static bool first_time = true;
    if (first_time) {
        first_time = false;
        std::ofstream str("C:/cygwin/home/stephen/projects/shallow-water/dump_to_file.txt");
    }

    const int START_FRAME = 0;  // inclusive (0=first frame, 1=second etc)
    const int STOP_FRAME = INT_MAX;   // exclusive
    static int frame_count = 0;
    if (frame_count >= START_FRAME && frame_count < STOP_FRAME) {

        std::ofstream str("C:/cygwin/home/stephen/projects/shallow-water/dump_to_file.txt", std::ios::app);

        str << "\n\nTIMESTEP NUMBER: " << frame_count << "\n";

        str << "\nInitial State (W, HU, HV, _):\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[sim_idx].get());
        
        str << "\nh (N/E/S/W):\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[1 - sim_idx].get());

        str << "\nu (N/E/S/W):\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[2].get());

        str << "\nv (N/E/S/W):\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[3].get());
    }
#endif

    // Pass 2
    // read: h, u, v
    // write: xflux, yflux

    vert_buf = m_psSimVertexBuffer10.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    ID3D11RenderTargetView * p2_tgt[] = { xflux_target, yflux_target, 0, 0 };
    context->OMSetRenderTargets(4, &p2_tgt[0], 0);

    context->PSSetShader(m_psSimPixelShader[1].get(), 0, 0);
    context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &u_tex);
    context->PSSetShaderResources(2, 1, &v_tex);

    context->Draw(6, 0);
    

#ifdef DUMP_TO_FILE
    if (frame_count >= START_FRAME && frame_count < STOP_FRAME) {
        std::ofstream str("C:/cygwin/home/stephen/projects/shallow-water/dump_to_file.txt", std::ios::app);

        str << "\nXFLUX: (W/HU/HV/_)\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[4].get());

        float wflux = GetSum(context, m_psFullSizeStagingTexture.get(), -1, -1, 0, ny-1);
        wflux -= GetSum(context, m_psFullSizeStagingTexture.get(), nx-1, nx-1, 0, ny-1);
        
        str << "\nYFLUX: (W/HU/HV/_)\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[5].get());

        wflux += GetSum(context, m_psFullSizeStagingTexture.get(), 0, nx-1, -1, -1);
        wflux -= GetSum(context, m_psFullSizeStagingTexture.get(), 0, nx-1, ny-1, ny-1);

        str << "\nBOUNDARY FLUX of w = " << wflux << std::endl;
    }
#endif

    
    // Pass 3
    // read: old_state, bottom, xflux, yflux
    // write: new_state

    vert_buf = m_psSimVertexBuffer00.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    context->PSSetShaderResources(0, 1, &pNULL); // unbind new_state_or_h_target so we can set it as output.

    ID3D11RenderTargetView * p3_tgt[] = { new_state_or_h_target, 0 };
    context->OMSetRenderTargets(2, &p3_tgt[0], 0);

    context->PSSetShader(m_psSimPixelShader[2].get(), 0, 0);
    context->PSSetShaderResources(0, 1, &old_state_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    context->PSSetShaderResources(2, 1, &xflux_tex);
    context->PSSetShaderResources(3, 1, &yflux_tex);

    context->Draw(6, 0);

#else

    // Lax Wendroff case

    // read: old state, bottom
    // write: new state

    vert_buf = m_psSimVertexBuffer00.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    ID3D11RenderTargetView * p3_tgt[] = { new_state_or_h_target, normal_target };
    context->OMSetRenderTargets(2, &p3_tgt[0], 0);

    context->PSSetShader(m_psLaxWendroffPixelShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &old_state_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);

    context->Draw(6, 0);

#endif
  
    
    // Boundary Conditions

    BoundaryConstBuffer cb;
    
    cb.boundary_epsilon = CalcEpsilon();
    cb.reflect_x = 2*nx+3;
    cb.reflect_y = 2*ny+3;
    cb.solid_wall_flag = (GetIntSetting("solid_walls") != 0);
    cb.sea_level = GetIntSetting("use_sea_level") ? GetSetting("sea_level") : -9999;

    const float inflow_width = 0.0;//GetSetting("inflow_width");
    const float inflow_height = 10.0;//GetSetting("inflow_height");
    const float W = GetSetting("valley_width");
    const float dx = GetSetting("valley_width") / (nx-1);
    cb.inflow_x_min = int((g_inlet_x - inflow_width + W/2)/dx) + 1;
    cb.inflow_x_max = int((g_inlet_x + inflow_width + W/2)/dx) + 3;
    cb.inflow_height = inflow_height;
    cb.inflow_speed = 0.01f;
    cb.g = GetSetting("gravity");
    cb.total_time = total_time;
    SeaWaveSettings('1', cb.sa1, cb.skx1, cb.sky1, cb.so1);
    SeaWaveSettings('2', cb.sa2, cb.skx2, cb.sky2, cb.so2);
    SeaWaveSettings('3', cb.sa3, cb.skx3, cb.sky3, cb.so3);
    SeaWaveSettings('4', cb.sa4, cb.skx4, cb.sky4, cb.so4);
    cb.sdecay = 0.01f / ny * GetSetting("valley_length");

    //memset(&cb, 0, sizeof(BoundaryConstBuffer));
    context->UpdateSubresource(m_psBoundaryConstantBuffer.get(), 0, 0, &cb, 0, 0);

    ID3D11Buffer *b_cst_buf = m_psBoundaryConstantBuffer.get();
    context->PSSetConstantBuffers(0, 1, &b_cst_buf);
    for (int i = 0; i < 3; ++i) context->PSSetShaderResources(1+i, 1, &pNULL);

    // use XFLUX as scratch space, then we'll copy back to the main output 
    context->OMSetRenderTargets(1, &xflux_target, 0);

    // now rebind the input as the output from the previous step (ie the new state).
    context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    
    for (int i = 0; i < 4; ++i) {
        vert_buf = m_psBoundaryVertexBuffer[i].get();
        context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

        context->PSSetShader(m_psBoundaryPixelShader[i].get(), 0, 0);
        
        context->Draw(6, 0);
    }
    
    // copy the temporary stuff from "xflux" back into the main texture.
    D3D11_BOX src_box;
    src_box.left = 2;
    src_box.right = nx + 2;
    src_box.top = 200;
    src_box.bottom = 250;
    src_box.front = 0;
    src_box.back = 1;
    context->CopySubresourceRegion(m_psSimTexture[1 - sim_idx].get(),
                                   0,  // subresource
                                   2,  // dest x
                                   200,   // dest y
                                   0,  // dest z
                                   m_psSimTexture[4].get(),  // xflux tex
                                   0, // subresource
                                   &src_box);

    src_box.top = ny+2;
    src_box.bottom = ny+4;
    context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, 2, ny+2, 0, m_psSimTexture[4].get(), 0, &src_box);

    src_box.left = 0;
    src_box.right = 2;
    src_box.top = 2;
    src_box.bottom = ny + 2;
    context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, 0, 2, 0, m_psSimTexture[4].get(), 0, &src_box);

    src_box.left = nx+2;
    src_box.right = nx+4;
    context->CopySubresourceRegion(m_psSimTexture[1-sim_idx].get(), 0, nx+2, 2, 0, m_psSimTexture[4].get(), 0, &src_box);


#ifdef DUMP_TO_FILE
    if (frame_count >= START_FRAME && frame_count < STOP_FRAME) {
        std::ofstream str("C:/cygwin/home/stephen/projects/shallow-water/dump_to_file.txt", std::ios::app);
        str << "\nFinal State:\n";
        DumpToFile(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[1-sim_idx].get());
        WriteTotal(str, context, m_psFullSizeStagingTexture.get(), m_psSimTexture[1-sim_idx].get());
    }
    ++frame_count;
#endif


#ifdef RUN_CHECKS
    RunChecks(context, m_psFullSizeStagingTexture.get(), m_psSimTexture[1-sim_idx].get());
#endif

#ifdef USE_KP07
    
    // Now do "pass 1" again, this means the H, U, V textures will be ready for the
    // next timestep. 
    // Also, the Normal texture will be created at this point.

    // first need to unbind 'old_state' from the pixel shader (as it is the target for our 'H' texture)
    context->PSSetShaderResources(0, 1, &pNULL);

    context->PSSetConstantBuffers(0, 1, &cst_buf);
    
    vert_buf = m_psSimVertexBuffer11.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    ID3D11RenderTargetView * p1_tgt[] = { old_state_or_h_target, u_target, v_target, normal_target };
    context->OMSetRenderTargets(4, &p1_tgt[0], 0);
    context->PSSetShader(m_psSimPixelShader[0].get(), 0, 0);
    context->PSSetShaderResources(0, 1, &new_state_or_h_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    context->Draw(6, 0);

#endif
    

    // Swap buffers.
    sim_idx = 1 - sim_idx;

    total_time += current_timestep;

    D3DPERF_EndEvent();

}

void ShallowWaterEngine::resetTimestep(float dt)
{
    D3DPERF_BeginEvent(0, L"demo_resettimestep");

    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
        
    // Run the GetStats pass
    // note: this uses the xflux, yflux textures as scratch space.

    ID3D11Buffer * cst_buf = m_psSimConstantBuffer.get();

    ID3D11ShaderResourceView * bottom_tex = m_psBottomTextureView.get();
    ID3D11ShaderResourceView * old_state_tex = m_psSimTextureView[sim_idx].get();
    ID3D11RenderTargetView * render_targets[] = { m_psSimRenderTargetView[4].get(),     // xflux texture
                                                  m_psSimRenderTargetView[5].get(),     // yflux texture
                                                  m_psGetStatsRenderTargetView.get() }; // dedicated R32_FLOAT texture
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psSimInputLayout.get());

    context->VSSetShader(m_psSimVertexShader.get(), 0, 0);

    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(nx + 4);
    vp.Height = float(ny + 4);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    context->VSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetConstantBuffers(0, 1, &cst_buf);

    ID3D11Buffer *vert_buf;
    const UINT stride = 8;
    const UINT offset = 0;

    vert_buf = m_psGetStatsVertexBuffer.get();
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    context->IASetIndexBuffer(m_psGetStatsIndexBuffer.get(), DXGI_FORMAT_R32_UINT, 0);
    
    context->OMSetRenderTargets(3, render_targets, 0);

    context->PSSetShader(m_psGetStatsPixelShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &old_state_tex);
    context->PSSetShaderResources(1, 1, &bottom_tex);
    
    //context->Draw(6, 0);
    context->DrawIndexed(6, 0, 0);

    // immediately read back the results.
    // (TODO: this will cause a pipeline stall... could run in a background thread to avoid this?)

    D3D11_BOX src_box;
    src_box.left = 0;
    src_box.right = nx/4;
    src_box.top = 0;
    src_box.bottom = ny/4;
    src_box.front = 0;
    src_box.back = 1;

    // copy first target into left part of the staging texture,
    // second target into right part.
    // third target goes into a separate staging texture.
    
    context->CopySubresourceRegion(m_psGetStatsStagingTexture4.get(),
                                   0,
                                   0,
                                   0,
                                   0,
                                   m_psSimTexture[4].get(),
                                   0,
                                   &src_box);

    context->CopySubresourceRegion(m_psGetStatsStagingTexture4.get(),
                                   0,
                                   nx/4,
                                   0,
                                   0,
                                   m_psSimTexture[5].get(),
                                   0,
                                   &src_box);

    context->CopySubresourceRegion(m_psGetStatsStagingTexture1.get(),
                                   0,
                                   0,
                                   0,
                                   0,
                                   m_psGetStatsTexture.get(),
                                   0,
                                   &src_box);
                          
    float mass = 0;
    float x_mtm = 0, y_mtm = 0;
    float ke = 0, pe = 0;
    float max_speed = 0, max_depth = 0;
    float cfl = 0;
    float max_froude = 0.0f;

    {
        MapTexture m(*context, *m_psGetStatsStagingTexture4);

        for (int j = 0; j < ny/4; ++j) {
            const char *row_ptr = reinterpret_cast<const char *>(m.msr.pData) + j * m.msr.RowPitch;
            for (int i = 0; i < nx/4; ++i) {
                const float *col_ptr_1 = reinterpret_cast<const float*>(row_ptr) + i * 4;
                const float *col_ptr_2 = reinterpret_cast<const float*>(row_ptr) + (i+nx/4) * 4;
                
                mass += col_ptr_1[0];   // sum(h)
                pe += col_ptr_1[1];     // sum(B*h + 0.5 * h^2)
                x_mtm += col_ptr_1[2];  // sum(hu)
                y_mtm += col_ptr_1[3];  // sum(hv)
                ke += col_ptr_2[0];     // sum(h*(u2+v2))
                max_speed = std::max(max_speed, col_ptr_2[1]);   // max(u2+v2)
                max_depth = std::max(max_depth, col_ptr_2[2]);   // max(h)
                cfl = std::max(cfl, col_ptr_2[3]);  // max((|u|+c)/dx, (|v|+c)/dy)
            }
        }
    }

    {
        MapTexture m(*context, *m_psGetStatsStagingTexture1);

        for (int j = 0; j < ny/4; ++j) {
            const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + j * m.msr.RowPitch;
            for (int i = 0; i < nx/4; ++i) {
                const float *p = reinterpret_cast<const float*>(row_ptr) + i;
                max_froude = std::max(max_froude, *p);   // max((u2+v2)/h)
            }
        }
    }
    
    const float DENSITY = 1000;   // kg m^-3
    const float AREA = GetSetting("valley_width") / float(nx-1) 
        * GetSetting("valley_length") / float(ny-1);   // m^2 (area of one cell)

    const float g = GetSetting("gravity");
    mass *= DENSITY * AREA;
    x_mtm *= DENSITY * AREA;
    y_mtm *= DENSITY * AREA;
    ke *= 0.5f * DENSITY * AREA;
    pe *= g * DENSITY * AREA;
    max_speed = std::sqrt(max_speed);
    max_froude /= g;
    max_froude = std::sqrt(max_froude);
        
    // The CFL number is cfl * dt, and this must be less than safety_factor, so dt < safety_factor/cfl
    const float safety_factor = GetSetting("max_cfl_number");
    current_timestep = std::min(dt * GetSetting("time_acceleration"), safety_factor / cfl);

    // update the displays
    SetSetting("mass", mass);
    SetSetting("x_momentum", x_mtm);
    SetSetting("y_momentum", y_mtm);
    SetSetting("kinetic_energy", ke);
    SetSetting("potential_energy", pe);
    SetSetting("total_energy", ke + pe);
    //SetSetting("testable", dt * GetSetting("time_acceleration"));
    //SetSetting("testable1", safety_factor / cfl);
    SetSetting("max_speed", max_speed);
    SetSetting("max_depth", max_depth);
    SetSetting("max_froude_number", max_froude);
    SetSetting("timestep", current_timestep);
    SetSetting("cfl_number", cfl * current_timestep);
    SetSetting("time_ratio", current_timestep / dt);
    
    fillConstantBuffers();  // communicate new dt to the simulation.

    D3DPERF_EndEvent();
}

void ShallowWaterEngine::render(ID3D11RenderTargetView *render_target_view)
{
    D3DPERF_BeginEvent(0, L"demo_render");

    // Setup rendering state
    context->ClearState();

    ID3D11ShaderResourceView * heightfield_tex = m_psTerrainTextureView.get();
    ID3D11ShaderResourceView * water_tex = m_psSimTextureView[sim_idx].get();   // {w, hu, hv, unused}
    ID3D11ShaderResourceView * normal_tex = m_psSimTextureView[6].get();     // {nX, nY, nZ, unused}
    ID3D11ShaderResourceView * skybox_tex = m_psSkyboxView.get();
    ID3D11ShaderResourceView * grass_tex = m_psGrassTextureView.get();
    ID3D11SamplerState * linear_sampler = m_psLinearSamplerState.get();
    
    
    // input assembler
    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psInputLayout.get());

    ID3D11Buffer *vert_buf = m_psMeshVertexBuffer.get();
    UINT stride = sizeof(MeshVertex);
    const UINT offset = 0;
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    context->IASetIndexBuffer(m_psMeshIndexBuffer.get(), DXGI_FORMAT_R32_UINT, 0);

    // vertex shader
    context->VSSetShader(m_psWaterVertexShader.get(), 0, 0);

    context->VSSetShaderResources(0, 1, &heightfield_tex);
    context->VSSetShaderResources(1, 1, &water_tex);
    context->VSSetShaderResources(2, 1, &normal_tex);
    
    ID3D11Buffer * cst_buf = m_psConstantBuffer.get();
    context->VSSetConstantBuffers(0, 1, &cst_buf);

    // viewport
    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(vp_width);
    vp.Height = float(vp_height);
    vp.MinDepth = 0;
    vp.MaxDepth = 1;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    context->RSSetViewports(1, &vp);

    // pixel shader
    context->PSSetShader(m_psWaterPixelShader.get(), 0, 0);
    context->PSSetShaderResources(0, 1, &grass_tex);
    context->PSSetShaderResources(1, 1, &skybox_tex);

    context->PSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetSamplers(0, 1, &linear_sampler);
    
    // setup depth buffer if we need to
    int rw, rh;
    int dw = -999, dh;
    GetRenderTargetSize(render_target_view, rw, rh);
    if (m_psDepthStencil.get()) GetTextureSize(m_psDepthStencil.get(), dw, dh);
    if (rw != dw || rh != dh) createDepthStencil(rw, rh);

    // clear the depth buffer
    context->ClearDepthStencilView(m_psDepthStencilView.get(), D3D11_CLEAR_DEPTH, 1.0f, 0);

    // set render target
    context->OMSetRenderTargets(1, &render_target_view, m_psDepthStencilView.get());

    // draw the mesh
    const int mesh_width = GetIntSetting("mesh_size_x");
    const int mesh_height = GetIntSetting("mesh_size_y");
    context->DrawIndexed(6 * (mesh_width - 1) * (mesh_height - 1), 0, 0);

    // now setup for the terrain mesh
    // (it uses the same mesh, but different shaders.)
    context->VSSetShader(m_psTerrainVertexShader.get(), 0, 0);
    context->PSSetShader(m_psTerrainPixelShader.get(), 0, 0);

    // draw the mesh
    context->DrawIndexed(6 * (mesh_width - 1) * (mesh_height - 1), 0, 0);


    // Now draw the skybox
    context->IASetInputLayout(m_psSkyboxInputLayout.get());
    vert_buf = m_psSkyboxVertexBuffer.get();
    stride = 3 * sizeof(float);   // 3 coords per vertex
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);
    context->VSSetShader(m_psSkyboxVertexShader.get(), 0, 0);
    context->PSSetShader(m_psSkyboxPixelShader.get(), 0, 0);
    context->Draw(30, 0);

    D3DPERF_EndEvent();
}


void ShallowWaterEngine::createShadersAndInputLayout()
{
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
    
    CreateVertexShader(device, "shallow_water.fx", "TerrainVertexShader", bytecode, m_psTerrainVertexShader);
    CreatePixelShader(device, "shallow_water.fx", "TerrainPixelShader", m_psTerrainPixelShader);

    CreateVertexShader(device, "shallow_water.fx", "WaterVertexShader", bytecode, m_psWaterVertexShader);
    CreatePixelShader(device, "shallow_water.fx", "WaterPixelShader", m_psWaterPixelShader);

    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "INXX_POSITION", 0, DXGI_FORMAT_R32G32_SINT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    ID3D11InputLayout *input_layout = 0;
    HRESULT hr = device->CreateInputLayout(&layout[0],
                                           1,
                                           bytecode->GetBufferPointer(),
                                           bytecode->GetBufferSize(),
                                           &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psInputLayout.reset(input_layout);


    CreateVertexShader(device, "kp07.hlsl", "SimVertexShader", bytecode, m_psSimVertexShader);
    CreatePixelShader(device, "kp07.hlsl", "Pass1", m_psSimPixelShader[0]);
    CreatePixelShader(device, "kp07.hlsl", "Pass2", m_psSimPixelShader[1]);
    CreatePixelShader(device, "kp07.hlsl", "Pass3", m_psSimPixelShader[2]);
    CreatePixelShader(device, "kp07.hlsl", "GetStats", m_psGetStatsPixelShader);
    CreatePixelShader(device, "kp07.hlsl", "NorthBoundary", m_psBoundaryPixelShader[0]);
    CreatePixelShader(device, "kp07.hlsl", "EastBoundary", m_psBoundaryPixelShader[1]);
    CreatePixelShader(device, "kp07.hlsl", "SouthBoundary", m_psBoundaryPixelShader[2]);
    CreatePixelShader(device, "kp07.hlsl", "WestBoundary", m_psBoundaryPixelShader[3]);

#ifndef USE_KP07
    CreatePixelShader(device, "lax_wendroff.hlsl", "LaxWendroffSinglePass", m_psLaxWendroffPixelShader);
#endif
    
    D3D11_INPUT_ELEMENT_DESC sim_layout[] = {
        { "TEX_IDX", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    hr = device->CreateInputLayout(&sim_layout[0],
                                   1,
                                   bytecode->GetBufferPointer(),
                                   bytecode->GetBufferSize(),
                                   &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psSimInputLayout.reset(input_layout);
}

// create the heightfield vertex and index buffers.
// (needs to be called again every time the mesh size changes.)
void ShallowWaterEngine::createMeshBuffers()
{
    const int width = GetIntSetting("mesh_size_x");
    const int height = GetIntSetting("mesh_size_y");

    // fill in the array of vertex positions
    // (these are texture indices -- add two to avoid the ghost zones)

    boost::scoped_array<MeshVertex> vertices(new MeshVertex[width*height]);
    
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            MeshVertex &v = vertices[j * width + i];
            v.x = i + 2;
            v.y = j + 2;
        }
    }

    // create data for the index buffer
    boost::scoped_array<int> indices(new int[6 * (width-1) * (height-1)]);
    int *p = &indices[0];
    for (int j = 0; j < height - 1; ++j) {
        for (int i = 0; i < width - 1; ++i) {
            const int tl = (j+1) * width + i;
            const int bl = j * width + i;
            const int tr = (j+1) * width + i+1;
            const int br = j * width + i+1;

            // write triangles in clockwise order.
            *p++ = tl;
            *p++ = br;
            *p++ = bl;

            *p++ = tl;
            *p++ = tr;
            *p++ = br;
        }
    }
            
    // create the vertex buffer
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(MeshVertex) * width * height;
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &vertices[0];

    ID3D11Buffer *pBuffer;
    HRESULT hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create mesh vertex buffer", hr);
    }
    m_psMeshVertexBuffer.reset(pBuffer);

    // create the index buffer
    bd.ByteWidth = sizeof(int) * 6 * (width-1) * (height-1);
    bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
    sd.pSysMem = &indices[0];
    sd.SysMemPitch = 0;

    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create mesh index buffer", hr);
    }
    m_psMeshIndexBuffer.reset(pBuffer);
}


void ShallowWaterEngine::createSimBuffers()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    
    CreateSimBuffer(device, m_psSimVertexBuffer11, 1, 1, nx + 3, ny + 3);   // single ghost layer around each side
    CreateSimBuffer(device, m_psSimVertexBuffer10, 1, 1, nx + 2, ny + 2);   // west/south ghost layer only
    CreateSimBuffer(device, m_psSimVertexBuffer00, 2, 2, nx + 2, ny + 2);   // interior zones only

    //CreateSimBuffer(device, m_psGetStatsVertexBuffer, 0, 0, nx/4, ny/4);
    CreateSimBuffer(device, m_psGetStatsVertexBuffer, m_psGetStatsIndexBuffer, 0, 0, nx/4, ny/4);

    CreateSimBuffer(device, m_psBoundaryVertexBuffer[0], 2, ny+2, nx+2, ny+4);  // north border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[1], nx+2, 2, nx+4, ny+2);  // east border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[2], 2, 200, nx+2, 250);   // south border
    CreateSimBuffer(device, m_psBoundaryVertexBuffer[3], 0, 2, 2, ny+2);   // west border
}

// Creates the terrain texture (leaving it uninitialized)
void ShallowWaterEngine::createTerrainTexture()
{
    // allow space for two ghost zones around each edge (four in total)
    CreateTexture(device,
                  GetIntSetting("mesh_size_x") + 4,
                  GetIntSetting("mesh_size_y") + 4,
                  0,  // initial data
                  DXGI_FORMAT_R32G32B32_FLOAT,
                  false,  // staging
                  m_psTerrainTexture,
                  &m_psTerrainTextureView,
                  0);

    CreateTexture(device,
                  GetIntSetting("mesh_size_x") + 4,
                  GetIntSetting("mesh_size_y") + 4,
                  0,  // initial data
                  DXGI_FORMAT_R32G32B32_FLOAT,
                  false, 
                  m_psBottomTexture,
                  &m_psBottomTextureView,
                  0);
}

// Creates and initializes the simulation textures
// Precondition: terrain heightfield is up to date
void ShallowWaterEngine::createSimTextures(ResetType reset_type)
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const float dam_pos = GetSetting("dam_position");

    const float xmin = -W/6;
    const float xmax = W/6;
    const float ymin = L/3;
    const float ymax = 2*L/3;

    /*
    float init_w = 0;
    if (reset_type == R_VALLEY) {
        // find the height of the lowest point along the dam
        float B_min = 99999999.f;
        for (float x = -W/2; x < W/2; x += W/(nx-1)) {
            B_min = std::min(B_min, GetTerrainHeight(x, dam_pos));
        }

        init_w = B_min + 1;
    } else if (reset_type == R_SEA) {
        init_w = GetSetting("sea_level");
    }
    */
            
    boost::scoped_array<float> ic(new float[(nx+4) * (ny+4) * 4]);
    float *p = &ic[0];

    for (int j = 0; j < ny + 4; ++j) {

        for (int i = 0; i < nx + 4; ++i) {

            const int ii = std::max(2, std::min(nx+1, i));
            const int jj = std::max(2, std::min(ny+1, j));
            const float B = g_bottom[(nx+4) * jj + ii].BA;

            const float x = (ii-2)*W/(nx-1) - W/2;
            const float y = (jj-2)*L/(ny-1);

            float w = B;
            /*
            if (reset_type == R_VALLEY) {
                if (y > dam_pos && B < init_w) {
                    w = init_w;

                    // add a wave pattern for some extra interest
                    w += 0.02f * sin(-0.2f*x + 0.8f*y);
                }
            } else if (reset_type == R_SEA) {
                w = std::max(B, init_w);
            } else if (reset_type == R_SQUARE && x >= xmin && x < xmax && y >= ymin && y < ymax) {
                w = B + 7.0f;
            }
            */

            // initial condition

            //@@ Temp
            if( j > (ny + 4)/4 &&
                j < (ny + 4)/2 && 
                i > 0 &&
                i < (nx + 4)/2 )
            {
                *p++ = (13.0);  // w
            }
            else
            {
                *p++ = w;  // w
            }
            *p++ = 0;  // hu
            *p++ = 0;  // hv
            *p++ = 0;  // unused
        }
    }
    
    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &ic[0];
    sd.SysMemPitch = (nx+4) * 4 * sizeof(float);

    for (int i = 0; i < 7; ++i) {
        CreateTexture(device,
                      nx + 4,
                      ny + 4,
                      i < 2 ? &sd : 0,
                      DXGI_FORMAT_R32G32B32A32_FLOAT,
                      false,
                      m_psSimTexture[i],
                      &m_psSimTextureView[i],
                      &m_psSimRenderTargetView[i]);
    }

    // TODO: This has no need to be (nx+4) by (ny+4),
    // we only ever use the top left quarter of it ((nx/4) by (ny/4))...
    // (Perhaps could do two or three separate drawcalls in GetStats pass instead of one multiple-output drawcall.)
    CreateTexture(device,
                  nx + 4,
                  ny + 4,
                  0,
                  DXGI_FORMAT_R32_FLOAT,
                  false,
                  m_psGetStatsTexture,
                  0,
                  &m_psGetStatsRenderTargetView);
    
    sim_idx = 0;
    bootstrap_needed = true;
    
    // Create a staging texture so we can read the data back again when required
    // (e.g. for debugging, or when changing terrain level)
    // TODO: we should be able to get away without using this staging texture; that would save a bit of memory
    CreateTexture(device,
                  nx + 4,
                  ny + 4,
                  0,
                  DXGI_FORMAT_R32G32B32A32_FLOAT,
                  true,
                  m_psFullSizeStagingTexture,
                  0,
                  0);

    // Create another staging texture of size (NX/2) * (NY/4) for GetStats
    CreateTexture(device,
                  nx/2,
                  ny/4,
                  0,
                  DXGI_FORMAT_R32G32B32A32_FLOAT,
                  true,
                  m_psGetStatsStagingTexture4,
                  0,
                  0);

    // and another of size (NX/4) * (NY/4) with only one channel
    CreateTexture(device,
                  nx/4,
                  ny/4,
                  0,
                  DXGI_FORMAT_R32_FLOAT,
                  true,
                  m_psGetStatsStagingTexture1,
                  0,
                  0);
}



// Updates the terrain texture given current settings.
// Also updates the water state texture such that the water depth remains unchanged.
// NOTE: If size has changed then call createTerrainTexture first.
void ShallowWaterEngine::fillTerrainTexture()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    // Get the existing water state onto the CPU
    context->CopyResource(m_psFullSizeStagingTexture.get(), m_psSimTexture[sim_idx].get());

    // Now loop through and change w values into h values
    // TODO: it would probably be better to do this on the GPU (would avoid a copy / copy back; and could
    // save memory for the staging texture as well).
    {
        MapTexture m(*context, *m_psFullSizeStagingTexture);
        for (int j = 0; j < ny+4; ++j) {
            char * row_ptr = reinterpret_cast<char*>(m.msr.pData) + j * m.msr.RowPitch;
            const BottomEntry *B_row_ptr = &g_bottom[j * (nx+4)];
            
            for (int i = 0; i < nx+4; ++i) {
                float *col_ptr = reinterpret_cast<float*>(row_ptr) + i * 4;
                const BottomEntry *B_col_ptr = B_row_ptr + i;
                
                (*col_ptr) -= (B_col_ptr->BA);
            }
        }
        
        // Compute the new terrain heightfield
        UpdateTerrainHeightfield();
        
        // Write the new terrain textures
        context->UpdateSubresource(m_psTerrainTexture.get(),
                                   0,   // subresource
                                   0,   // overwrite whole resource
                                   &g_terrain_heightfield[0],
                                   (nx+4) * 12,
                                   0);  // slab pitch
        context->UpdateSubresource(m_psBottomTexture.get(),
                                   0,  // subresource
                                   0,  // overwrite whole resource
                                   &g_bottom[0],
                                   (nx+4) * 12,
                                   0); // slab pitch
        
        // Loop through and change h values back into w values
        for (int j = 0; j < ny+4; ++j) {
            char * row_ptr = reinterpret_cast<char*>(m.msr.pData) + j * m.msr.RowPitch;
            const BottomEntry *B_row_ptr = &g_bottom[j * (nx+4)];
            
            for (int i = 0; i < nx+4; ++i) {
                float *col_ptr = reinterpret_cast<float*>(row_ptr) + i * 4;
                const BottomEntry *B_col_ptr = B_row_ptr + i;
                
                (*col_ptr) += ( (B_col_ptr->BA) );
            }
        }
    }
    
    // Copy water texture back to the GPU
    context->CopyResource(m_psSimTexture[sim_idx].get(), m_psFullSizeStagingTexture.get());

    // need to re-bootstrap
    bootstrap_needed = true;
}



// Updates the terrain texture, but does not touch the water texture
void ShallowWaterEngine::fillTerrainTextureLite()
{
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    // Compute the new terrain heightfield
    UpdateTerrainHeightfield();
 
    std::ofstream stream(L"H:/Water/TestData/dump_to_file_shallow.txt", std::ios::app);
    stream << "\n\n mesh_size_x: " << nx << "\n";
    stream << "mesh_size_y: " << ny << "\n";
    stream << "g_terrain_heightfield: (intx, inty, B, dBdx, dBdy)" << "\n";
    for(int i=0;i<ny+4;i++)
    {
        for(int j=0;j<nx+4;j++)
        {
            stream << "(" << j << ", " << i << ", " << g_terrain_heightfield[i*nx+j].B << ", " << g_terrain_heightfield[i*nx+j].dBdx << ", " << g_terrain_heightfield[i*nx+j].dBdy << ") \n";
        }

    }
    stream << "\n\n g_bottom: (intx, inty, BA, BX, BY)" << "\n";
    for(int i=0;i<ny+4;i++)
    {
        for(int j=0;j<nx+4;j++)
        {
            stream << "(" << j << ", " << i << ", " << g_bottom[i*nx+j].BA << ", " << g_bottom[i*nx+j].BX << ", " << g_bottom[i*nx+j].BY << ") \n";
        }

    }

    // Write the new terrain textures
    context->UpdateSubresource(m_psTerrainTexture.get(),
                               0,   // subresource
                               0,   // overwrite whole resource
                               &g_terrain_heightfield[0],
                               (nx+4) * 12,
                               0);  // slab pitch
    context->UpdateSubresource(m_psBottomTexture.get(),
                               0,  // subresource
                               0,  // overwrite whole resource
                               &g_bottom[0],
                               (nx+4) * 12,
                               0); // slab pitch
    
    // need to re-bootstrap
    bootstrap_needed = true;
}



// create the constant buffers, leaving them uninitialized
void ShallowWaterEngine::createConstantBuffers()
{
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = RoundUpTo16(sizeof(MyConstBuffer));
    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;

    ID3D11Buffer *pBuffer;
    HRESULT hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psConstantBuffer.reset(pBuffer);

    bd.ByteWidth = RoundUpTo16(sizeof(SimConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psSimConstantBuffer.reset(pBuffer);

    bd.ByteWidth = RoundUpTo16(sizeof(BoundaryConstBuffer));
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psBoundaryConstantBuffer.reset(pBuffer);
}

// update the constant buffers given current settings
void ShallowWaterEngine::fillConstantBuffers()
{
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    MyConstBuffer cb;

    const float sun_alt = GetSetting("sun_alt") * PI / 180.0f;
    const float sun_az = GetSetting("sun_az") * PI / 180.0f;
    cb.light_dir = XMFLOAT3(cos(sun_alt)*sin(sun_az), cos(sun_alt)*cos(sun_az), sin(sun_alt));
    cb.ambient = GetSetting("ambient");

    cb.eye_mult = XMFLOAT3(-W/(nx-1), -L/(ny-1), -1.0f);
    cb.eye_trans = XMFLOAT3(W/2 + 2*W/(nx-1) + camera_x, 2*L/(ny-1) + camera_y, camera_z);

    cb.world_mult_x = W/(nx-1);
    cb.world_mult_y = L/(ny-1);
    cb.world_trans_x = -W/2 - 2*W/(nx-1);
    cb.world_trans_y = -2*L/(ny-1);

    const float grass_size = 15.0f; // width of the entire BMP in metres
    cb.terrain_tex_scale = XMFLOAT2(W / (nx-1) / grass_size, L / (ny-1) / grass_size);
    cb.world_to_grass_tex_x = 1.0f / grass_size;
    cb.world_to_grass_tex_y = 1.0f / grass_size;
    cb.grass_tex_of_origin_x = W/2 / grass_size;
    cb.grass_tex_of_origin_y = 0;

    cb.fresnel_coeff = GetSetting("fresnel_coeff");
    cb.fresnel_exponent = GetSetting("fresnel_exponent");
    cb.specular_intensity = GetSetting("specular_intensity");
    cb.specular_exponent = GetSetting("specular_exponent");
    cb.refractive_index = GetSetting("refractive_index");
    cb.attenuation_1 = GetSetting("attenuation_1");
    cb.attenuation_2 = GetSetting("attenuation_2");
    cb.deep_r = GetSetting("deep_r");
    cb.deep_g = GetSetting("deep_g");
    cb.deep_b = GetSetting("deep_b");
    cb.nx_plus_1 = nx + 1;
    cb.ny_plus_1 = ny + 1;

    // setup matrix transforming (tex_x, tex_y, world_z, 1) into normalized device coordinates.

    //const XMMATRIX tex_to_world( W / (nx - 1), 0,             0,  -W/2 - 2*W/(nx-1),
    //                             0,            L / (ny - 1),  0,       - 2*L/(ny-1),
    //                             0,            0,             1,  0,
    //                             0,            0,             0,  1  );
    const XMMATRIX scale_world( W / (nx - 1), 0,             0,  0,
                                 0,            L / (ny - 1),  0,  0,
                                 0,            0,             1,  0,
                                 0,            0,             0,  1  );

    const XMMATRIX translate_world( 1,            0,             0,         -W/2 - 2*W/(nx-1),
                                    0,            1,             0,         - 2*L/(ny-1),
                                    0,            0,             1,         0,
                                    0,            0,             0,         1  );
    
    const XMMATRIX translate_camera_pos( 1, 0, 0, -camera_x,
                                         0, 1, 0, -camera_y,
                                         0, 0, 1, -camera_z,
                                         0, 0, 0, 1 );

    const XMMATRIX rotate_around_z( cos(yaw), -sin(yaw), 0, 0,
                                    sin(yaw), cos(yaw),  0, 0,
                                    0,        0,         1, 0,
                                    0,        0,         0, 1 );

    const XMMATRIX rotate_around_x( 1, 0,           0,          0,
                                    0, cos(pitch),  sin(pitch), 0,
                                    0, -sin(pitch), cos(pitch), 0,
                                    0, 0,           0,          1 );

    const XMMATRIX swap_y_and_z( 1, 0, 0, 0,
                                 0, 0, 1, 0,
                                 0, 1, 0, 0,
                                 0, 0, 0, 1 );

    const float far = far_plane_dist;
    const float near = near_plane_dist;
    const XMMATRIX perspective( xpersp, 0,      0,                  0,
                                0,      ypersp, 0,                  0,
                                0,      0,      far / (far - near), -near * far / (far - near),
                                0,      0,      1,                  0  );

    //cb.tex_to_clip = perspective * swap_y_and_z * rotate_around_x
    //    * rotate_around_z * translate_camera_pos * tex_to_world;
    cb.tex_to_clip = perspective * swap_y_and_z * rotate_around_x
        * rotate_around_z * translate_camera_pos * translate_world * scale_world;

    cb.skybox_mtx = perspective * swap_y_and_z * rotate_around_x * rotate_around_z;

    // Now write it to the const buffer
    context->UpdateSubresource(m_psConstantBuffer.get(),
                               0,     // subresource
                               0,     // overwrite whole resource
                               &cb,
                               0,     // row pitch
                               0);    // slab pitch


    // Now do the sim parameters
    SimConstBuffer sb;

    const float theta = GetSetting("theta");  // 1=most dissipative, 2=most oscillatory, 1.3 = good default.
    const float g = GetSetting("gravity");
    const float dt = current_timestep;
    
    const float dx = W / (nx-1);
    const float dy = L / (ny-1);
    
    sb.two_theta = 2 * theta;
    sb.two_over_nx_plus_four = 2.0f / (nx+4);
    sb.two_over_ny_plus_four = 2.0f / (ny+4);
    sb.g = g;
    sb.half_g = 0.5f * g;
    sb.g_over_dx = g / dx;
    sb.g_over_dy = g / dy;
    sb.one_over_dx = 1.0f / dx;
    sb.one_over_dy = 1.0f / dy;
    sb.dt = dt;
    sb.epsilon = CalcEpsilon();
    sb.nx = nx;
    sb.ny = ny;
    sb.friction = GetSetting("friction");

    // Now write it to the constant buffer
    context->UpdateSubresource(m_psSimConstantBuffer.get(),
                               0,  // subresource
                               0,  // overwrite whole buffer
                               &sb,
                               0,  // row pitch
                               0); // slab pitch
}

void ShallowWaterEngine::createDepthStencil(int w, int h)
{
    D3D11_TEXTURE2D_DESC td;
    memset(&td, 0, sizeof(td));
    td.Width = w;
    td.Height = h;
    td.MipLevels = 1;
    td.ArraySize = 1;
    td.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    td.SampleDesc.Count = 1;
    td.Usage = D3D11_USAGE_DEFAULT;
    td.BindFlags = D3D11_BIND_DEPTH_STENCIL;

    ID3D11Texture2D *pDepthStencil = 0;
    HRESULT hr = device->CreateTexture2D(&td, 0, &pDepthStencil);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create depth stencil texture", hr);
    }
    m_psDepthStencil.reset(pDepthStencil);
    
    D3D11_DEPTH_STENCIL_VIEW_DESC vd;
    memset(&vd, 0, sizeof(vd));
    vd.Format = td.Format;
    vd.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;

    ID3D11DepthStencilView *pDepthStencilView = 0;
    hr = device->CreateDepthStencilView(pDepthStencil, &vd, &pDepthStencilView);
    if (FAILED(hr)) {
        throw Coercri::DXError("failed to create depth stencil view", hr);
    }
    m_psDepthStencilView.reset(pDepthStencilView);
}

void ShallowWaterEngine::loadGraphics()
{
    std::ifstream str("grass01.bmp", std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    
    CreateTextureWithMips(device,
                          parr->getWidth(),
                          parr->getHeight(),
                          reinterpret_cast<unsigned char *>(&(*parr)(0,0)),
                          m_psGrassTexture,
                          m_psGrassTextureView);

    D3D11_SAMPLER_DESC sd;
    memset(&sd, 0, sizeof(sd));
    sd.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sd.AddressU = sd.AddressV = sd.AddressW = D3D11_TEXTURE_ADDRESS_WRAP;
    sd.ComparisonFunc = D3D11_COMPARISON_NEVER;
    sd.MaxLOD = D3D11_FLOAT32_MAX;

    ID3D11SamplerState * sampler = 0;
    HRESULT hr = device->CreateSamplerState(&sd, &sampler);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create sampler state", hr);
    }
    m_psLinearSamplerState.reset(sampler);
}

/*
void ShallowWaterEngine::createBlendState()
{
    D3D11_BLEND_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.RenderTarget[0].BlendEnable = true;
    bd.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
    bd.RenderTarget[0].DestBlend = D3D11_BLEND_SRC_ALPHA;
    bd.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
    bd.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ZERO;
    bd.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ZERO;
    bd.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
    bd.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

    ID3D11BlendState *blend_state;
    HRESULT hr = device->CreateBlendState(&bd, &blend_state);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create blend state", hr);
    }
    m_psWaterBlendState.reset(blend_state);
}
*/

void ShallowWaterEngine::loadSkybox()
{
    // load the cube texture
    ID3D11ShaderResourceView * srv;
    HRESULT hr = D3DX11CreateShaderResourceViewFromFile(device, "skybox_c.dds", 0, 0, &srv, 0);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to load the skybox", hr);
    }
    m_psSkyboxView.reset(srv);

    // create shaders
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
    CreateVertexShader(device, "shallow_water.fx", "SkyboxVertexShader", bytecode, m_psSkyboxVertexShader);
    CreatePixelShader(device, "shallow_water.fx", "SkyboxPixelShader", m_psSkyboxPixelShader);

    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };

    // create vertex buffer input layout
    ID3D11InputLayout *input_layout = 0;
    hr = device->CreateInputLayout(&layout[0],
                                   1,
                                   bytecode->GetBufferPointer(),
                                   bytecode->GetBufferSize(),
                                   &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout (skybox) failed", hr);
    }
    m_psSkyboxInputLayout.reset(input_layout);

    // create vertex buffer

    const float vertices[5*6*3] = {

        // north (+Y)
        -0.5f, 0.5f, 0.5f,
        0.5f,  0.5f, 0.5f,
        0.5f,  0.5f, -0.5f,
        0.5f,  0.5f, -0.5f,
        -0.5f, 0.5f, -0.5f,
        -0.5f, 0.5f, 0.5f,

        // west (-X)
        -0.5f, -0.5f, 0.5f,
        -0.5f, 0.5f,  0.5f,
        -0.5f, 0.5f,  -0.5f,
        -0.5f, 0.5f,  -0.5f,
        -0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, 0.5f,

        // south (-Y)
        0.5f, -0.5f, 0.5f,
        -0.5f, -0.5f, 0.5f,
        -0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, 0.5f,

        // east (+X)
        0.5f, 0.5f, 0.5f,
        0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, 0.5f, -0.5f,
        0.5f, 0.5f, 0.5f,
        
        // up (+Z)
        -0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f, 0.5f,
        -0.5f, -0.5f, 0.5f,
    };
   
    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(vertices);
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &vertices[0];

    ID3D11Buffer *pBuffer;
    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create vertex buffer for skybox", hr);
    }
    m_psSkyboxVertexBuffer.reset(pBuffer);
}


// mouse picking

namespace {
    const int MOUSE_PICK_NUM_TRIANGLES = 12;

    struct LeftMouseVertex {
        float x, y;
    };
}

// converts pixel position into a world point on the water surface

// returns true if there was a hit; in this case world_xyz will contain the hit point.
// otherwise, returns false (and corrupts world_xyz).

bool ShallowWaterEngine::mousePick(int mouse_x, int mouse_y,
                                   float &world_x, float &world_y, float &world_z, float &depth, float &u, float &v)
{
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");

    const XMMATRIX translate_camera_pos( 1, 0, 0, camera_x,
                                         0, 1, 0, camera_y,
                                         0, 0, 1, camera_z,
                                         0, 0, 0, 1 );
    const XMMATRIX rotate_around_z( cos(yaw), sin(yaw), 0, 0,
                                    -sin(yaw), cos(yaw), 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 1 );
    const XMMATRIX rotate_around_x( 1, 0, 0, 0,
                                    0, cos(pitch), -sin(pitch), 0,
                                    0, sin(pitch), cos(pitch), 0,
                                    0, 0, 0, 1);
    const XMMATRIX swap_y_and_z( 1, 0, 0, 0,
                                 0, 0, 1, 0,
                                 0, 1, 0, 0,
                                 0, 0, 0, 1 );

    // transpose because D3D treats position vectors as row vectors for some reason...
    const XMMATRIX eye_to_world = XMMatrixTranspose(translate_camera_pos * rotate_around_z * 
        rotate_around_x * swap_y_and_z);

    const float ndc_x = 2 * (mouse_x + 0.5f) / vp_width - 1;
    const float ndc_y = -2 * (mouse_y + 0.5f) / vp_height + 1;

    const float stepsize = 1;  // in metres

    bool found = false;
    XMFLOAT4 world_pos;


    // step along the ray in multiples of stepsize.
    // this is a bit inefficient, but it works well enough for what we want to use it for...

    // this assumes resetTimestep has previously been called.
    MapTexture m(*context, *m_psGetStatsStagingTexture4);

    for (float eye_z = 1; eye_z < 1000 && !found; eye_z += stepsize) {

        XMFLOAT4 eye_pos( eye_z * ndc_x / xpersp,
                          eye_z * ndc_y / ypersp,
                          eye_z,
                          1 );
        XMVECTOR eye_pos_vec = XMLoadFloat4(&eye_pos);
        XMVECTOR world_pos_vec = XMVector4Transform(eye_pos_vec, eye_to_world);
        XMStoreFloat4(&world_pos, world_pos_vec);
        
        // clip to the terrain region
        if (world_pos.x >= -W/2 && world_pos.x <= W/2 
        && world_pos.y >= 0 && world_pos.y <= L) {

            // convert world pos to a texture position in the staging texture.
            const float tex_x = (world_pos.x + W/2) / W;
            const float tex_y = (world_pos.y) / L;
            int ix = int(tex_x * (nx/4));
            int iy = int(tex_y * (ny/4));
            if (ix<0) ix=0;
            if (iy<0) iy=0;
            if (ix>nx/4-1) ix=nx/4-1;
            if (iy>ny/4-1) iy=ny/4-1;

            // lookup the "h" value at this point
            const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + iy * m.msr.RowPitch;
            const float *col_ptr = reinterpret_cast<const float*>(row_ptr) + ix * 4;
            const float avg_h = (*col_ptr) / 16.0f;

            // lookup terrain height at this position
            const float B = GetTerrainHeight(world_pos.x, world_pos.y);
            const float w = B + avg_h;
            
            if (world_pos.z <= w) {
                world_x = world_pos.x;
                world_y = world_pos.y;
                world_z = w;
                depth = avg_h;

                u = CalcU(avg_h, col_ptr[2]/16);
                v = CalcU(avg_h, col_ptr[3]/16);
                
                found = true;
            }
        }
    }

    return found;
}

void ShallowWaterEngine::applyMouseShader(float world_x, float world_y, float dt)
{
    const int action = GetIntSetting("left_mouse_action");
    if (action == LM_RAISE_TERRAIN || action == LM_LOWER_TERRAIN) {
        raiseLowerTerrain(world_x, world_y, dt);
        return;
    }
    
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const float brush_radius = GetSetting("left_mouse_radius");
    const float strength = GetSetting("left_mouse_strength") * dt;
    
    // transformation from world coords to texture index.
    const float world_to_idx_x_scale = (nx-1)/W;
    const float world_to_idx_x_bias = (nx+4)/2.0f;
    const float world_to_idx_y_scale = (ny-1)/L;
    const float world_to_idx_y_bias = 2.5f;

    // transformation from "radius 1" coords (in the vertex buf) to texture indices.
    const float x_scale = world_to_idx_x_scale * brush_radius;
    const float x_bias = world_to_idx_x_scale * world_x + world_to_idx_x_bias;
    const float y_scale = world_to_idx_y_scale * brush_radius;
    const float y_bias = world_to_idx_y_scale * world_y + world_to_idx_y_bias;

    // setup the const buffer
    LeftMouseConstBuffer cb;
    cb.scale_x = x_scale;
    cb.scale_y = y_scale;
    cb.bias_x = x_bias;
    cb.bias_y = y_bias;
    cb.two_over_nx_plus_four = 2.0f / (nx+4);
    cb.two_over_ny_plus_four = 2.0f / (ny+4);

    switch (action) {
    case LM_ADD_WATER:
        cb.disp_A = strength;
        cb.disp_B = 0;
        break;

    case LM_REMOVE_WATER:
        cb.disp_A = -strength;
        cb.disp_B = 0;
        break;

    case LM_STIR_WATER:
        cb.disp_A = -strength;
        cb.disp_B = 1.5f * strength;
        break;

    default:
        cb.disp_A = cb.disp_B = 0;
        break;
    }
    
    context->UpdateSubresource(m_psLeftMouseConstantBuffer.get(),
                               0, 0, &cb, 0, 0);

    // prepare to copy state texture to a temporary work area.
    // find the texture indices for the copy region.
    const int centre_ix = int(world_to_idx_x_scale * world_x + world_to_idx_x_bias);
    const int centre_iy = int(world_to_idx_y_scale * world_y + world_to_idx_y_bias);
    const int radius_x = int(world_to_idx_x_scale * brush_radius) + 2;  // safety margin
    const int radius_y = int(world_to_idx_y_scale * brush_radius) + 2;

    const int ix_min = std::max(0, std::min(nx+4, centre_ix - radius_x));   // inclusive
    const int ix_max = std::max(0, std::min(nx+4, centre_ix + radius_x + 1));  // exclusive
    const int iy_min = std::max(0, std::min(ny+4, centre_iy - radius_y));
    const int iy_max = std::max(0, std::min(ny+4, centre_iy + radius_y + 1));
    
    // do the copy
    D3D11_BOX src_box;
    src_box.left = ix_min;
    src_box.right = ix_max;
    src_box.top = iy_min;
    src_box.bottom = iy_max;
    src_box.front = 0;
    src_box.back = 1;
    context->CopySubresourceRegion(m_psSimTexture[4].get(),  // dest texture (xflux -- being used here as scratch space)
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psSimTexture[sim_idx].get(),  // src texture
                                   0,  // subresource
                                   &src_box);

    // now we can set up the pixel shader to input from xflux texture,
    // and write back to current state texture.
    
    context->ClearState();

    context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->IASetInputLayout(m_psLeftMouseInputLayout.get());

    ID3D11Buffer *vert_buf = m_psLeftMouseVertexBuffer.get();
    const UINT stride = 8;
    const UINT offset = 0;
    context->IASetVertexBuffers(0, 1, &vert_buf, &stride, &offset);

    ID3D11Buffer *cst_buf = m_psLeftMouseConstantBuffer.get();
    context->VSSetShader(m_psLeftMouseVertexShader.get(), 0, 0);
    context->VSSetConstantBuffers(0, 1, &cst_buf);

    D3D11_VIEWPORT vp;
    memset(&vp, 0, sizeof(vp));
    vp.Width = float(nx+4);
    vp.Height = float(ny+4);
    vp.MaxDepth = 1;
    context->RSSetViewports(1, &vp);

    ID3D11ShaderResourceView *tex_views[] = { m_psSimTextureView[4].get(), m_psBottomTextureView.get() };
    context->PSSetShader(m_psLeftMousePixelShader.get(), 0, 0);
    context->PSSetConstantBuffers(0, 1, &cst_buf);
    context->PSSetShaderResources(0, 2, &tex_views[0]);

    ID3D11RenderTargetView * rtv = m_psSimRenderTargetView[sim_idx].get();
    context->OMSetRenderTargets(1, &rtv, 0);

    context->Draw(MOUSE_PICK_NUM_TRIANGLES * 3, 0);
}

void ShallowWaterEngine::raiseLowerTerrain(float world_x, float world_y, float dt)
{
    // Do this on CPU as it is the easiest way -- just want to get this code written
    // quickly, don't care about efficiency at this point :)

    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    const float dx = W / (nx-1);
    const float dy = L / (ny-1);    
    const float brush_radius = GetSetting("left_mouse_radius") * 2; // *2 to account for 'gaussian' nature of brush
    float strength = GetSetting("left_mouse_strength");
    const int action = GetIntSetting("left_mouse_action");    
    float displacement;
    if (action == LM_LOWER_TERRAIN) displacement = -strength*dt; else displacement = strength*dt;
    
    const float world_to_idx_x_scale = (nx-1)/W;
    const float world_to_idx_x_bias = (nx+4)/2.0f;
    const float world_to_idx_y_scale = (ny-1)/L;
    const float world_to_idx_y_bias = 2.5f;

    const int centre_ix = int(world_to_idx_x_scale * world_x + world_to_idx_x_bias);
    const int centre_iy = int(world_to_idx_y_scale * world_y + world_to_idx_y_bias);
    const int radius_x = int(world_to_idx_x_scale * brush_radius) + 3;  // safety margin
    const int radius_y = int(world_to_idx_y_scale * brush_radius) + 3;

    const int ix_min = std::max(0, std::min(nx+4, centre_ix - radius_x));   // inclusive
    const int ix_max = std::max(0, std::min(nx+4, centre_ix + radius_x + 1));  // exclusive
    const int iy_min = std::max(0, std::min(ny+4, centre_iy - radius_y));
    const int iy_max = std::max(0, std::min(ny+4, centre_iy + radius_y + 1));    

    D3D11_BOX src_box;
    src_box.left = ix_min;
    src_box.right = ix_max;
    src_box.top = iy_min;
    src_box.bottom = iy_max;
    src_box.front = 0;
    src_box.back = 1;
    context->CopySubresourceRegion(m_psFullSizeStagingTexture.get(),
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psSimTexture[sim_idx].get(),  // src texture -- current state
                                   0,  // subresource
                                   &src_box);

    {
        MapTexture m(*context, *m_psFullSizeStagingTexture);

        for (int iy = iy_min; iy < iy_max; ++iy) {

            // have to be careful about the exact definitions of BX, BY, BA (see terrain_heightfield.cpp)
            // (Water level w must be brought up in line with BA)
            
            const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
            const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;

            char * row_ptr = static_cast<char*>(m.msr.pData) + m.msr.RowPitch * iy;
            
            for (int ix = ix_min; ix < ix_max; ++ix) {

                const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
                const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;

                const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
                const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
                const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
                const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;

                float * ptr = reinterpret_cast<float*>(row_ptr) + 4 * ix;

                *ptr += (plus_plus + plus_minus + minus_plus + minus_minus) / 4;
            }
        }
    }

    context->CopySubresourceRegion(m_psSimTexture[sim_idx].get(),  // dest texture
                                   0,  // subresource
                                   ix_min,  // dest x
                                   iy_min,  // dest y
                                   0,  // dest z
                                   m_psFullSizeStagingTexture.get(),  // src texture
                                   0,  // subresource
                                   &src_box);    

    // Now we are going to update the in-memory copy of BottomTexture and
    // upload it to the GPU via UpdateSubresource.

    for (int iy = iy_min; iy < iy_max; ++iy) {

        const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;

        BottomEntry * row_ptr = &g_bottom[iy * (nx+4)];
            
        for (int ix = ix_min; ix < ix_max; ++ix) {
            
            const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            
            const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
            const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
            const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
            const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;

            BottomEntry * ptr = row_ptr + ix;

            ptr->BX += (plus_plus + plus_minus)/2;
            ptr->BY += (plus_plus + minus_plus)/2;
            ptr->BA += (plus_plus + plus_minus + minus_plus + minus_minus)/4;
        }
    }

    context->UpdateSubresource(m_psBottomTexture.get(),  // dest texture
                               0,  // subresource
                               &src_box,  // dest box
                               &g_bottom[iy_min * (nx+4) + ix_min],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused)

    // more crappy copy paste code...
    
    for (int iy = iy_min; iy < iy_max; ++iy) {

        const float wy_plus = (iy + 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        const float wy_minus = (iy - 0.5f - world_to_idx_y_bias) / world_to_idx_y_scale - world_y;
        
        TerrainEntry * row_ptr = &g_terrain_heightfield[iy * (nx+4)];
            
        for (int ix = ix_min; ix < ix_max; ++ix) {

            const float wx_plus = (ix + 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;
            const float wx_minus = (ix - 0.5f - world_to_idx_x_bias) / world_to_idx_x_scale - world_x;

            const float plus_plus = GaussianBrush(wx_plus, wy_plus, brush_radius) * displacement;
            const float plus_minus = GaussianBrush(wx_plus, wy_minus, brush_radius) * displacement;
            const float minus_plus = GaussianBrush(wx_minus, wy_plus, brush_radius) * displacement;
            const float minus_minus = GaussianBrush(wx_minus, wy_minus, brush_radius) * displacement;
            
            TerrainEntry * ptr = row_ptr + ix;

            ptr->B += (plus_plus + plus_minus + minus_plus + minus_minus)/4;
            ptr->dBdx += (plus_plus + plus_minus - minus_plus - minus_minus) / (2*dx);
            ptr->dBdy += (plus_plus + minus_plus - plus_minus - minus_minus) / (2*dy);
        }
    }  

    context->UpdateSubresource(m_psTerrainTexture.get(),  // dest texture
                               0,  // subresource
                               &src_box,  // dest box
                               &g_terrain_heightfield[iy_min * (nx+4) + ix_min],  // src data
                               (nx+4) * 3 * sizeof(float),   // row pitch
                               0);  // depth pitch (unused)    
}

void ShallowWaterEngine::setupMousePicking()
{
    // create shaders
    Coercri::ComPtrWrapper<ID3DBlob> bytecode;
    CreateVertexShader(device, "left_mouse.hlsl", "LeftMouseVertexShader", bytecode, m_psLeftMouseVertexShader);
    CreatePixelShader(device, "left_mouse.hlsl", "LeftMousePixelShader", m_psLeftMousePixelShader);
        
    // create input layout
    D3D11_INPUT_ELEMENT_DESC layout[] = {
        { "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
    };
    
    ID3D11InputLayout *input_layout = 0;
    HRESULT hr = device->CreateInputLayout(&layout[0],
                                           1,
                                           bytecode->GetBufferPointer(),
                                           bytecode->GetBufferSize(),
                                           &input_layout);
    if (FAILED(hr)) {
        throw Coercri::DXError("CreateInputLayout failed", hr);
    }
    m_psLeftMouseInputLayout.reset(input_layout);

    // create a vertex buffer for a circle of radius 1.
    std::vector<LeftMouseVertex> buf(MOUSE_PICK_NUM_TRIANGLES * 3);
    const float d_theta = 2*PI/MOUSE_PICK_NUM_TRIANGLES;
    for (int i = 0; i < MOUSE_PICK_NUM_TRIANGLES; ++i) {
        const float theta = d_theta * i;
        const float theta2 = d_theta * (i+1);

        buf[3*i].x = 0;
        buf[3*i].y = 0;

        buf[3*i+1].x = cos(theta);
        buf[3*i+1].y = sin(theta);

        buf[3*i+2].x = cos(theta2);
        buf[3*i+2].y = sin(theta2);
    }

    D3D11_BUFFER_DESC bd;
    memset(&bd, 0, sizeof(bd));
    bd.ByteWidth = sizeof(LeftMouseVertex) * MOUSE_PICK_NUM_TRIANGLES * 3;
    bd.Usage = D3D11_USAGE_IMMUTABLE;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;

    D3D11_SUBRESOURCE_DATA sd;
    memset(&sd, 0, sizeof(sd));
    sd.pSysMem = &buf[0];

    ID3D11Buffer *pBuffer = 0;
    hr = device->CreateBuffer(&bd, &sd, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create left mouse vertex buffer", hr);
    }
    m_psLeftMouseVertexBuffer.reset(pBuffer);

    // create the constant buffer.
    bd.ByteWidth = RoundUpTo16(sizeof(LeftMouseConstBuffer));
    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    hr = device->CreateBuffer(&bd, 0, &pBuffer);
    if (FAILED(hr)) {
        throw Coercri::DXError("Failed to create constant buffer", hr);
    }
    m_psLeftMouseConstantBuffer.reset(pBuffer);
}

float ShallowWaterEngine::getWaterHeight(float world_x, float world_y)
{
    /*
    const float W = GetSetting("valley_width");
    const float L = GetSetting("valley_length");
    const int nx = GetIntSetting("mesh_size_x");
    const int ny = GetIntSetting("mesh_size_y");
    MapTexture m(*context, *m_psGetStatsStagingTexture4);
    const float tex_x = (world_x + W/2) / W;
    const float tex_y = (world_y) / L;
    int ix = int(tex_x * (nx/4));
    int iy = int(tex_y * (ny/4));
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    const char *row_ptr = reinterpret_cast<const char*>(m.msr.pData) + iy * m.msr.RowPitch;
    const float *col_ptr = reinterpret_cast<const float*>(row_ptr) + ix * 4;
    const float avg_h = (*col_ptr) / 16.0f;
    */
    const float B = GetTerrainHeight(world_x, world_y);
    return 0.0;
    //return B + avg_h;
}
