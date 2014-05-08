/*
 * FILE:
 *   main.cpp
 *
 * PURPOSE:
 *   Main function for the Shallow Water demo
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

#include "engine.hpp"
#include "gui_manager.hpp"
#include "settings.hpp"
#include "terrain_heightfield.hpp"

#include "coercri/dx11/gfx/dx11_gfx_driver.hpp"
#include "coercri/dx11/gfx/dx11_window.hpp"
#include "coercri/gfx/gfx_context.hpp"
#include "coercri/gfx/window_listener.hpp"
#include "coercri/timer/generic_timer.hpp"

#include "boost/scoped_ptr.hpp"

#include <cmath>

const int GUI_WIDTH = 400;
int g_width = 500, g_height = 500;
bool g_quit = false;
bool g_resize = false;

namespace {
    const float PI = std::atan(1.0f) * 4.0f;

    void UpdateMousePick(int mx, int my, ShallowWaterEngine &engine)
    {
        float wx, wy, wz, d, u, v;
        if (mx >= 0 && mx < g_width && my >= 0 && my < g_height
        && engine.mousePick(mx, my, wx, wy, wz, d, u, v)) {
            SetSetting("x", wx);
            SetSetting("y", wy);
            SetSetting("z", wz);
            SetSetting("depth", d);
            SetSetting("u", u);
            SetSetting("v", v);

            const float spd = std::sqrt(u*u + v*v);
            const float c = std::sqrt(GetSetting("gravity") * d);
            SetSetting("froude", spd / c);
            
        } else {
            SetSetting("x", -1);
            SetSetting("y", -1);
            SetSetting("z", -1);
            SetSetting("depth", -1);
            SetSetting("u", -1);
            SetSetting("v", -1);
            SetSetting("froude", -1);
        }
    }
}

class MyListener : public Coercri::WindowListener {
public:
    MyListener(ShallowWaterEngine &e, GuiManager &g, Coercri::Window &w) : 
        engine(e), gui_manager(g), window(w),
        fwd(false), bwd(false), left(false), right(false), up(false), down(false),
        shift(false),
        cam_x(0), cam_y(0), cam_z(30),
        pitch(0), yaw(0),
        mx(-999), my(-999), left_mouse_down(false), right_mouse_down(false)
    { }

    void onLoseFocus() { fwd = bwd = left = right = up = down = shift = false; left_mouse_down = right_mouse_down = false; }
    
    void onClose() { g_quit = true; }
    void onResize(int w, int h) {
        const int gw = gui_manager.isGuiShown() ? GUI_WIDTH : 0;
        g_width = std::max(0, w - gw);
        g_height = h;
        g_resize = true;
    }

    void onRawKey(bool pressed, Coercri::RawKey rk)
    {
        switch (rk) {
        case Coercri::RK_W: fwd = pressed; break;
        case Coercri::RK_S: bwd = pressed; break;
        case Coercri::RK_A: left = pressed; break;
        case Coercri::RK_D: right = pressed; break;
        case Coercri::RK_PAGE_UP: up = pressed; break;
        case Coercri::RK_PAGE_DOWN: down = pressed; break;
        case Coercri::RK_LEFT_SHIFT: case Coercri::RK_RIGHT_SHIFT: shift = pressed; break;
        }
    }

    void onMouseDown(int x, int y, Coercri::MouseButton mb)
    {
        if (mb == Coercri::MB_RIGHT) {
            mx = x;
            my = y;
            right_mouse_down = true;
            window.captureMouse(true);

        } else if (mb == Coercri::MB_LEFT) {
            left_mouse_down = true;
        }
    }

    void onMouseUp(int x, int y, Coercri::MouseButton mb)
    {
        if (mb == Coercri::MB_RIGHT) {
            right_mouse_down = false;
            window.captureMouse(false);
        } else if (mb == Coercri::MB_LEFT) {
            left_mouse_down = false;
        }
    }

    void onMouseMove(int new_x, int new_y) 
    {
        if (right_mouse_down) {
            const int dx = new_x - mx;
            const int dy = new_y - my;

            const float angle_speed = 0.001f;

            yaw += dx * angle_speed;
            pitch -= dy * angle_speed;

            if (pitch > PI/2) pitch = PI/2;
            if (pitch < -PI/2) pitch = -PI/2;
        }

        mx = new_x;
        my = new_y;
    }

    void update(float dt)
    {
        // HACK to communicate camera resets back from the gui.
        gui_manager.getCameraReset(cam_x, cam_y, cam_z, pitch, yaw);
        
        // move the camera
        float speed = shift ? 300.0f : 30.0f;
        speed *= dt;

        const float W = GetSetting("valley_width");
        const float L = GetSetting("valley_length");

        if (fwd) {
            cam_x += sin(yaw) * cos(pitch) * speed;
            cam_y += cos(yaw) * cos(pitch) * speed;
            cam_z += sin(pitch) * speed;
        }
        if (bwd) {
            cam_x -= sin(yaw) * speed;
            cam_y -= cos(yaw) * speed;
            cam_z -= sin(pitch) * speed;
        }
        if (left) {
            cam_x -= cos(yaw) * speed;
            cam_y += sin(yaw) * speed;
        }
        if (right) {
            cam_x += cos(yaw) * speed;
            cam_y -= sin(yaw) * speed;
        }
        if (up) {
            cam_z += cos(pitch) * speed;
            cam_x -= sin(pitch) * sin(yaw) * speed;
            cam_y -= sin(pitch) * cos(yaw) * speed;
        }
        if (down) {
            cam_z -= cos(pitch) * speed;
            cam_x += sin(pitch) * sin(yaw) * speed;
            cam_y += sin(pitch) * cos(yaw) * speed;
        }

        // clip camera position
        const float max_ratio = 0.5f;
        float max_clip = max_ratio*std::max(W,L);
        if (GetIntSetting("clip_camera")) max_clip = 0;
        if (cam_x < -W/2 -max_clip) cam_x = -W/2 - max_clip;
        if (cam_x > W/2 + max_clip) cam_x = W/2 + max_clip;
        if (cam_y < -max_clip) cam_y = -max_clip;
        if (cam_y > L + max_clip) cam_y = L + max_clip;

        const float height = engine.getWaterHeight(std::min(std::max(-W/2, cam_x), W/2),
                                                   std::min(std::max(0.0f, cam_y), L));
        cam_z = std::max(1.5f + height, cam_z);
        cam_z = std::min(500.0f + height, cam_z);        
        
        // TODO: only need to call this if the camera actually moved
        // (or the window resized...)
        engine.moveCamera(cam_x, cam_y, cam_z, yaw, pitch, g_width, g_height);

        // do left mouse interaction
        if (left_mouse_down) {
            const float wy = GetSetting("y");
            if (wy >= 0) {
                const float wx = GetSetting("x");
                const float wy = GetSetting("y");
                engine.applyMouseShader(wx, wy, dt);
            }

            //// rainfall - constant intensity
            //if( wy >= 0 )
            //{
            //    const float width = GetSetting("valley_width");
            //    const float length = GetSetting("valley_length");
            //    const float cell = 20.0;
            //    const int countW = static_cast<int>( width / cell );
            //    const int countL = static_cast<int>( length / cell );
            //    const float StartW = -0.5 * width;
            //    const float StartL = 0.0;
            //    for(int i=0;i<countL;++i)
            //    {
            //        for(int j=0;j<countW;++j)
            //        {
            //            const float offsetCellW = static_cast<float>( rand() % (int)cell ); // within one cell
            //            const float offsetCellL = static_cast<float>( rand() % (int)cell ); // within one cell
            //            const float wx = StartW + j * cell + offsetCellW;
            //            const float wy = StartL + i * cell + offsetCellL;
            //            engine.applyMouseShader(wx, wy, dt);
            //        }
            //    }
            //    
            //}
        }
    }

    int getMX() const { return mx; }
    int getMY() const { return my; }

private:

    ShallowWaterEngine &engine;
    GuiManager &gui_manager;
    Coercri::Window &window;

    bool fwd, bwd, left, right, up, down;
    bool shift;
        
    float cam_x, cam_y, cam_z;
    float pitch, yaw;

    int mx, my;
    bool left_mouse_down, right_mouse_down;
};

int real_main()
{
    // Setup Direct3D
    
#ifdef _DEBUG
    const bool debug = true;
#else
    const bool debug = false;
#endif
    
    // using texture lookups in the vertex shader, therefore must have feature level >= 10.0.

    std::auto_ptr<Coercri::DX11GfxDriver> gfx_driver(new Coercri::DX11GfxDriver(D3D_DRIVER_TYPE_HARDWARE,
                                                                                D3D11_CREATE_DEVICE_SINGLETHREADED | (debug ? D3D11_CREATE_DEVICE_DEBUG : 0),
                                                                                D3D_FEATURE_LEVEL_10_0));
    boost::shared_ptr<Coercri::Timer> timer(new Coercri::GenericTimer);

    // set up the gui
    std::auto_ptr<MyListener> listener;
    boost::shared_ptr<Coercri::DX11Window> window = 
        boost::static_pointer_cast<Coercri::DX11Window>(
            gfx_driver->createWindow(g_width + GUI_WIDTH, g_height, true, false, "Shallow Water Demo - Copyright (C) Stephen Thompson 2012"));
    GuiManager gui_manager(window, timer, GUI_WIDTH);

    // Create the ShallowWaterEngine
    boost::scoped_ptr<ShallowWaterEngine> engine(
        new ShallowWaterEngine(gfx_driver->getDevice(), gfx_driver->getDeviceContext()));

    listener.reset(new MyListener(*engine, gui_manager, *window));
    window->addWindowListener(listener.get());
    

    unsigned int last_draw = timer->getMsec();
    bool is_gui_shown = true;
    
    g_resize = true; // make sure it "resizes" first time

    while (!g_quit) {

        // slight hack to make sure window resizes when gui is shown/hidden
        bool new_is_gui_shown = gui_manager.isGuiShown();
        if (new_is_gui_shown != is_gui_shown) {
            int win_w, win_h;
            window->getSize(win_w, win_h);
            listener->onResize(win_w, win_h);
            is_gui_shown = new_is_gui_shown;
        }
        
        if (g_resize) {
            gui_manager.resize();
            g_resize = false;
        }

        // see if updates need to be done
        switch (g_reset_type) {
        case R_MESH:
        case R_VALLEY:
        case R_SEA:
        case R_SQUARE:
            engine->remesh(g_reset_type);
            break;

        case R_TERRAIN:
            engine->newTerrainSettings();
            break;
        };

        g_reset_type = R_NONE;        
        
        engine->resetTimestep(0.001f);
        int timestep_count = 0;
        unsigned int timer_at_zero_timesteps = timer->getMsec();
        
        //temp
        const int steps_between_reset = 10;

        while (!g_resize && g_reset_type == R_NONE && !g_quit && gui_manager.isGuiShown() == is_gui_shown) 
        {
            //engine->resetTimestep(0.001f);

            const unsigned int frame_time = 1;  // in msec. acts as fps limiter.

            // see if dt needs to be reset
            const int steps_between_reset = 10;
            if (timestep_count >= steps_between_reset) {
                unsigned int timer_now = timer->getMsec();
                float dt = float(timer_now - timer_at_zero_timesteps) / (1000.0f * steps_between_reset);
                SetSetting("testable", timer_now - timer_at_zero_timesteps);
                SetSetting("testable1", dt);
                engine->resetTimestep(dt);
                timestep_count = 0;
                timer_at_zero_timesteps = timer_now;

                UpdateMousePick(listener->getMX(), listener->getMY(), *engine);
            }
            
            gfx_driver->pollEvents();
            gui_manager.logic();

            if (g_reset_type != R_NONE) break;  // don't continue if the settings are out of date
            
            const unsigned int time_now = timer->getMsec();
            const int time_since_last = int(time_now - last_draw);
            
            bool do_render = true; // int(time_now - last_draw) >= frame_time;
            if (window->needsRepaint()) do_render = true;
            
            if (do_render) {

                last_draw = time_now;

                // do update
                listener->update(float(time_since_last) / 1000.0f);

                // do timestep
                const int timesteps_per_frame = GetIntSetting("timesteps_per_frame");
                for (int i = 0; i < timesteps_per_frame; ++i) 
                {
                    engine->timestep();
                    ++timestep_count;
                }

                // clear the screen
                const float rgba[] = { 0, 0, 0, 1 };
                gfx_driver->getDeviceContext()->ClearRenderTargetView(window->getRenderTargetView(), rgba);

                // draw the landscape & water
                engine->render(window->getRenderTargetView());

                // draw the gui (using coercri 2D routines)
                {
                    std::auto_ptr<Coercri::GfxContext> gc = window->createGfxContext();
                    gui_manager.draw(*gc);
                }

                window->cancelInvalidRegion();
            
            } else {
                const int time_to_next_frame = int(last_draw + frame_time - time_now);
                timer->sleepMsec(std::max(1, time_to_next_frame + 1));
            }
        }
    }

    return 0;
}

int main()
{
    try {
        return real_main();
    } catch (std::exception &e) {
        MessageBox(0, e.what(), "Error", MB_ICONEXCLAMATION | MB_OK);
    } catch (...) {
        MessageBox(0, "Unknown exception", "Error", MB_ICONEXCLAMATION | MB_OK);
    }

    return 1;
}

int CALLBACK WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
    return main();
}
