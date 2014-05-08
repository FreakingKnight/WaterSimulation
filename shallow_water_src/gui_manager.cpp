/*
 * FILE:
 *   gui_manager.cpp
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
#include "gui_manager.hpp"

#include "coercri/gfx/bitmap_font.hpp"
#include "coercri/gfx/load_bmp.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

namespace {
    // Constantly creating and destroying ostringstream objects seems v. bad
    // for performance, so we create a global one and re-use it
    std::ostringstream g_stringstream;
    void SetLabelTxt(gcn::Label &label, double value)
    {
        g_stringstream.str("");
        g_stringstream << value;
        label.setCaption(g_stringstream.str());
        label.adjustSize();
    }
    
    double RealToSlider(const Setting &setting, double value)
    {
        /*
        if (setting.logarithmic) {
            return std::log(value);
        } else {
            return value;
        }
        */
        return value;
    }

    double SliderToReal(const Setting &setting, double value)
    {
        /*
        if (setting.logarithmic) {
            return std::exp(value);
        } else {
            return value;
        }
        */
        return value;
    }

    class MyListModel : public gcn::ListModel {
    public:
        std::string getElementAt(int i) {
            if (i >= 0 && i < int(elts.size())) return elts[i];
            else return "";
        }

        int getNumberOfElements() {
            return int(elts.size());
        }

        void add(const std::string &s) { elts.push_back(s); }

    private:
        std::vector<std::string> elts;
    };



    void SetupValley()
    {
        SetSettingD("mesh_size_x", 300);
        SetSettingD("mesh_size_y", 900);
        SetSettingD("solid_walls", 0);
        SetSettingD("inflow_width", 4);
        SetSettingD("inflow_height", 1);
        SetSettingD("gravity", 10);
        SetSettingD("friction", 0.02);
        SetSettingD("theta", 1.1);
        SetSettingD("max_cfl_number", 0.2);
        SetSettingD("timesteps_per_frame", 4);
        SetSettingD("time_acceleration", 1);
        SetSettingD("clip_camera", 1);
        SetSettingD("valley_length", 300);
        SetSettingD("valley_width", 100);
        SetSettingD("valley_wall_height", 6);
        SetSettingD("gradient_top", 0.1);
        SetSettingD("gradient_bottom", 0.02);
        SetSettingD("valley_shape", 2);
        SetSettingD("channel_depth_top", 2);
        SetSettingD("channel_depth_bottom", 2);
        SetSettingD("channel_width_top", 10);
        SetSettingD("channel_width_bottom", 10);
        SetSettingD("dam_on", 1);
        SetSettingD("dam_height", 10);
        SetSettingD("dam_position", 150);
        SetSettingD("dam_middle_width", 4);
        SetSettingD("dam_middle_height", -5);
        SetSettingD("dam_thickness", 8);
        SetSettingD("meander_wavelength", 100);
        SetSettingD("meander_amplitude", 8);
        SetSettingD("meander_fractal", 0.3);
        SetSettingD("use_sea_level", 0);
        SetSettingD("sea_level", 10);
        SetSettingD("sa1", 0.3);
        SetSettingD("sk1", 0.224);
        SetSettingD("sk1_dir", 5.18);
        SetSettingD("so1", 1);
        SetSettingD("sa2", 0.25);
        SetSettingD("sk2", 0.201);
        SetSettingD("sk2_dir", 4.83);
        SetSettingD("so2", 1.07);
        SetSettingD("sa3", 0.1);
        SetSettingD("sk3", 0.3);
        SetSettingD("sk3_dir", 4.9);
        SetSettingD("so3", 4);
        SetSettingD("sa4", 0);
        SetSettingD("sk4", 0);
        SetSettingD("sk4_dir", 0);
        SetSettingD("s04", 0);
        SetSettingD("fov", 45);
        SetSettingD("sun_alt", 25);
        SetSettingD("sun_az", 300);
        SetSettingD("ambient", 0.75);
        SetSettingD("fresnel_coeff", 0.983);
        SetSettingD("fresnel_exponent", 5);
        SetSettingD("specular_intensity", 1);
        SetSettingD("specular_exponent", 15);
        SetSettingD("refractive_index", 1.33);
        SetSettingD("attenuation_1", 0.08);
        SetSettingD("attenuation_2", 0.08);
        SetSettingD("deep_r", 0.05);
        SetSettingD("deep_g", 0.1);
        SetSettingD("deep_b", 0.2);
    }

    void SetupValleyHires()
    {
        SetupValley();
        SetSettingD("mesh_size_x", 600);
        SetSettingD("mesh_size_y", 1200);
        SetSettingD("inflow_width", 2);
        SetSettingD("timesteps_per_frame", 10);
        SetSettingD("valley_length", 100);
        SetSettingD("valley_width", 50);
        SetSettingD("channel_width_top", 5);
        SetSettingD("channel_width_bottom", 5);
        SetSettingD("dam_position", 50);
        SetSettingD("dam_height", 7);
        SetSettingD("dam_middle_height", -3);
        SetSettingD("dam_thickness", 4);
    }

    void SetupSea()
    {
        SetupValley();
        SetSettingD("mesh_size_x", 600);
        SetSettingD("mesh_size_y", 300);
        SetSettingD("inflow_width", 0);
        SetSettingD("inflow_height", 0);
        SetSettingD("valley_length", 300);
        SetSettingD("valley_width", 600);
        SetSettingD("valley_wall_height", 0);
        SetSettingD("channel_depth_top", 0);
        SetSettingD("channel_depth_bottom", 0);
        SetSettingD("channel_width_top", 0);
        SetSettingD("channel_width_bottom", 0);
        SetSettingD("dam_on", 0);
        SetSettingD("use_sea_level", 1);
        SetSettingD("sea_level", 8);
        SetSettingD("attenuation_2", 0.5);
        SetSettingD("deep_r", 0.01);
        SetSettingD("deep_g", 0.08);
        SetSettingD("deep_b", 0.1);
    }

    void SetupFlatPlane()
    {
        SetupValley();
        SetSettingD("mesh_size_x", 400);
        SetSettingD("mesh_size_y", 400);
        SetSettingD("solid_walls", 1);
        SetSettingD("inflow_width", 0);
        SetSettingD("inflow_height", 0);
        SetSettingD("valley_length", 200);
        SetSettingD("valley_width", 200);
        SetSettingD("valley_wall_height", 0);
        SetSettingD("gradient_top", 0);
        SetSettingD("gradient_bottom", 0);
        SetSettingD("channel_depth_top", 0);
        SetSettingD("channel_depth_bottom", 0);
        SetSettingD("channel_width_top", 0);
        SetSettingD("channel_width_bottom", 0);
        SetSettingD("dam_on", 0);
        SetSettingD("clip_camera", 0);
    }
}

GuiManager::GuiManager(boost::shared_ptr<Coercri::Window> window_,
                       boost::shared_ptr<Coercri::Timer> timer_,
                       int gui_width_)
    : gui(new gcn::Gui),
      window(window_),
      listener(window_, gui, timer_),
      gui_width(gui_width_),
      last_time(0), frame_count(-999),
      timer(timer_),
      gui_shown(true),
      cam_reset(false), cam_x(0), cam_y(0), cam_z(0), cam_pitch(0), cam_yaw(0)
{
    // load a font
    std::ifstream str("knights_sfont.bmp", std::ios::in | std::ios::binary);
    boost::shared_ptr<Coercri::PixelArray> parr = Coercri::LoadBMP(str);
    font.reset(new Coercri::BitmapFont(parr));
    cg_font.reset(new Coercri::CGFont(font));
    gcn::Widget::setGlobalFont(cg_font.get());

    // add the CGListener
    window->addWindowListener(&listener);
    listener.enableGui();

    // Create the widgets.
    createGui();

    // Default settings
    // note: duplicated code (see action())
    SetupValley();
    g_reset_type = R_VALLEY;
    cam_reset = true;
    cam_x = -90;
    cam_y = 280;
    cam_z = 40;
    cam_pitch = -0.2f;
    cam_yaw = 2.7f;
    resetSliders();    
}

GuiManager::~GuiManager()
{
    window->rmWindowListener(&listener);
}

void GuiManager::createGui()
{
    const int value_label_width = 80;

    show_gui_container.reset(new gcn::Container);
    hide_gui_button.reset(new gcn::Button("Hide GUI"));
    show_gui_button.reset(new gcn::Button("Show GUI"));
    hide_gui_button->addActionListener(this);
    show_gui_button->addActionListener(this);
    show_gui_container->add(show_gui_button.get());

    reset_simulation_button.reset(new gcn::Button("Reset Simulation"));
    reset_simulation_button->addActionListener(this);
    
    container.reset(new gcn::Container);

    tabbed_area.reset(new gcn::TabbedArea);
    container->add(tabbed_area.get());
    tabbed_area->setSize(container->getWidth(), container->getHeight());

    int max_label_width = 0;

    boost::shared_ptr<gcn::Container> presets_tab(new gcn::Container);
    presets_tab->setSize(container->getWidth(), container->getHeight());
    tabs.push_back(presets_tab);
    tabbed_area->addTab("Presets", presets_tab.get());

    const int pad = 14;
    int y = 12;
    preset_button_valley.reset(new gcn::Button("Valley"));
    presets_tab->add(preset_button_valley.get(), 5, y);
    preset_button_valley->addActionListener(this);
    y += pad + preset_button_valley->getHeight();
    
    preset_button_valley_hires.reset(new gcn::Button("Valley (high res)"));
    presets_tab->add(preset_button_valley_hires.get(), 5, y);
    preset_button_valley_hires->addActionListener(this);
    y += pad + preset_button_valley->getHeight();

    preset_button_flat.reset(new gcn::Button("Flat plane"));
    presets_tab->add(preset_button_flat.get(), 5, y);
    preset_button_flat->addActionListener(this);
    y += pad + preset_button_valley->getHeight();

    preset_button_sea.reset(new gcn::Button("Sea + waves"));
    presets_tab->add(preset_button_sea.get(), 5, y);
    preset_button_sea->addActionListener(this);

    const int pmaxw = std::max(std::max(std::max(preset_button_valley->getWidth(),
                                                 preset_button_valley_hires->getWidth()),
                                        preset_button_sea->getWidth()),
                               preset_button_flat->getWidth());
    preset_button_valley->setWidth(pmaxw);
    preset_button_valley_hires->setWidth(pmaxw);
    preset_button_sea->setWidth(pmaxw);
    preset_button_flat->setWidth(pmaxw);
    preset_button_valley->setAlignment(gcn::Graphics::LEFT);
    preset_button_valley_hires->setAlignment(gcn::Graphics::LEFT);
    preset_button_sea->setAlignment(gcn::Graphics::LEFT);
    preset_button_flat->setAlignment(gcn::Graphics::LEFT);
    
    
    for (const Setting *setting = &g_settings[0]; setting->name; ++setting) {

        if (setting->type == S_NEW_TAB && setting->name[0] != 0) {
            boost::shared_ptr<gcn::Container> current_tab(new gcn::Container);
            current_tab->setSize(container->getWidth(), container->getHeight());
            tabs.push_back(current_tab);
            tabbed_area->addTab(setting->name, current_tab.get());
        }
        
        boost::shared_ptr<gcn::Label> lab(new gcn::Label(setting->name));
        labels.push_back(lab);

        boost::shared_ptr<gcn::Slider> slider;
        boost::shared_ptr<gcn::CheckBox> checkbox;
        boost::shared_ptr<gcn::DropDown> dropdown;

        switch (setting->type) {
        case S_SLIDER:
        case S_SLIDER_MULT_4:
        case S_SLIDER_INT:
            slider.reset(new gcn::Slider(RealToSlider(*setting, setting->min),
                                         RealToSlider(*setting, setting->max)));
            slider->setValue(RealToSlider(*setting, setting->value));
            slider->addActionListener(this);
            break;

        case S_CHECKBOX:
            checkbox.reset(new gcn::CheckBox("", setting->value == 1));
            checkbox->addActionListener(this);
            break;

        case S_LEFT_MOUSE_DROPDOWN:

            {
                boost::shared_ptr<MyListModel> list_model(new MyListModel);
                list_model->add("Add Water");
                list_model->add("Remove Water");
                list_model->add("Stir Water");
                list_model->add("Raise Terrain");
                list_model->add("Lower Terrain");

                list_models.push_back(list_model);
            
                dropdown.reset(new gcn::DropDown(list_model.get()));
                dropdown->addActionListener(this);
            }
            break;
        }

        sliders.push_back(slider);
        check_boxes.push_back(checkbox);
        dropdowns.push_back(dropdown);

        boost::shared_ptr<gcn::Label> val_lbl;
        if (setting->type == S_LABEL || setting->type == S_SLIDER || setting->type == S_SLIDER_MULT_4 || setting->type == S_SLIDER_INT) {
            val_lbl.reset(new gcn::Label);
            SetLabelTxt(*val_lbl, setting->value);
        }
        values.push_back(val_lbl);
        
        const int w = lab->getWidth();
        if (w > max_label_width) max_label_width = w;
    }

    const int slider_width = std::max(10, gui_width - max_label_width - value_label_width);

    y = 0;

    std::vector<boost::shared_ptr<gcn::Container> >::iterator tab_it = tabs.begin();
    max_y = 0;
    
    for (size_t i = 0; i < sliders.size(); ++i) {

        if (g_settings[i].type == S_NEW_TAB) {
            ++tab_it;
            max_y = std::max(max_y, y);
            
            y = 8;
            if (tab_it == tabs.end()) {
                max_y += tabbed_area->getFont()->getHeight() + 30;
                y = max_y + 20;
            }
            
            continue;
        }
        
        if (sliders[i]) {
            sliders[i]->setWidth(slider_width);
            sliders[i]->setHeight(labels[i]->getHeight());
        }

        if (values[i]) {
            values[i]->setWidth(value_label_width);
            values[i]->setHeight(labels[i]->getHeight());
        }

        if (dropdowns[i]) {
            dropdowns[i]->setWidth(slider_width + value_label_width/2);
        }

        gcn::Container & con =
            (tab_it == tabs.end() ? *container : **tab_it);
        
        con.add(labels[i].get(), 0, y);
        if (sliders[i]) con.add(sliders[i].get(), max_label_width, y);
        if (check_boxes[i]) con.add(check_boxes[i].get(), max_label_width, y);
        if (dropdowns[i]) con.add(dropdowns[i].get(), max_label_width, y);

        int valx = max_label_width;
        if (g_settings[i].type != S_LABEL) valx += slider_width;
        if (values[i]) con.add(values[i].get(), valx, y);
        
        y += std::max(labels[i]->getHeight(),
                      dropdowns[i] ? dropdowns[i]->getHeight() : 0);
    }

    y += 16;
    container->add(reset_simulation_button.get(), 8, y);
    container->add(hide_gui_button.get(), 8 + reset_simulation_button->getWidth() + 50, y);
}

void GuiManager::resize()
{
    // gui_width is fixed, so we don't have to resize any of the widgets;
    // we just move the container to a new position.

    // (we also shrink the container if the window is so small that
    // gui_width doesn't fit.)
    
    int win_width, win_height;
    window->getSize(win_width, win_height);
    const int width = std::min(win_width, gui_width);
    const int x = win_width - width;
    
    container->setWidth(width);
    container->setHeight(win_height);
    container->setX(x);
    show_gui_container->setSize(win_width, win_height);
    show_gui_container->setPosition(0, 0);
    show_gui_container->setOpaque(false);

    tabbed_area->setSize(container->getWidth(), max_y);
    for (std::vector<boost::shared_ptr<gcn::Container> >::iterator it = tabs.begin(); it != tabs.end(); ++it) {
        (*it)->setSize(container->getWidth(), max_y);
    }

    show_gui_button->setPosition(win_width - show_gui_button->getWidth() - 10,
                                 10);
    show_gui_button->adjustSize();

    if (gui_shown) {
        gui->setTop(container.get());
    } else {
        gui->setTop(show_gui_container.get());
    }
}

void GuiManager::logic()
{
    listener.processInput();

    // reset all value labels (every few frames)
    if (frame_count == 0) {
        for (int i = 0; i < int(sliders.size()); ++i) {
            if (values[i]) {
                SetLabelTxt(*values[i], g_settings[i].value);
            }
        }
    }
}

void GuiManager::draw(Coercri::GfxContext &gc)
{
    listener.draw(gc);

    if (frame_count < 0) {
        frame_count = 0;
        last_time = timer->getMsec();

    } else {
        ++frame_count;
        unsigned int time_now = timer->getMsec();
        unsigned int difference = time_now - last_time;
        if (difference > 250) {

            const float fps = (1000.0f * frame_count) / float(difference);
            frame_count = 0;
            last_time = time_now;

            SetSetting("fps", fps);
        }
    }
}

void GuiManager::action(const gcn::ActionEvent &e)
{
    gcn::Slider *slid = static_cast<gcn::Slider*>(e.getSource());
    gcn::CheckBox *cb = static_cast<gcn::CheckBox*>(e.getSource());
    gcn::DropDown *dd = static_cast<gcn::DropDown*>(e.getSource());

    for (size_t i = 0; i < sliders.size(); ++i) {
        if (sliders[i].get() == slid) {
            g_settings[i].value = SliderToReal(g_settings[i], slid->getValue());
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);

            if (g_settings[i].type == S_SLIDER_MULT_4) {
                int val = int(g_settings[i].value + 0.5);
                val = (val + 3) & (~3);
                g_settings[i].value = double(val);
            } else if (g_settings[i].type == S_SLIDER_INT) {
                int val = int(g_settings[i].value + 0.5);
                g_settings[i].value = double(val);
            }
            
            SetLabelTxt(*values[i], g_settings[i].value);
            break;
        }

        if (check_boxes[i].get() == cb) {
            g_settings[i].value = cb->isSelected() ? 1.0 : 0.0;
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);
            break;
        }

        if (dropdowns[i].get() == dd) {
            g_settings[i].value = dd->getSelected();
            g_reset_type = std::max(g_reset_type, g_settings[i].reset_type);
            break;
        }
    }

    if (e.getSource() == hide_gui_button.get()) {
        gui_shown = false;
        resize();
    } else if (e.getSource() == show_gui_button.get()) {
        gui_shown = true;
        resize();
    } else if (e.getSource() == reset_simulation_button.get()) {
        g_reset_type = R_MESH;
    } else if (e.getSource() == preset_button_valley.get()) {
        // note: code duplicated in ctor
        SetupValley();
        g_reset_type = R_VALLEY;
        cam_reset = true;
        cam_x = -90;
        cam_y = 280;
        cam_z = 40;
        cam_pitch = -0.2f;
        cam_yaw = 2.7f;
        resetSliders();
    } else if (e.getSource() == preset_button_valley_hires.get()) {
        SetupValleyHires();
        g_reset_type = R_VALLEY;
        cam_reset = true;
        cam_x = -7;
        cam_y = 20;
        cam_z = 7;
        cam_pitch = -0.3f;
        cam_yaw = 0.0f;
        resetSliders();
    } else if (e.getSource() == preset_button_sea.get()) {
        SetupSea();
        g_reset_type = R_SEA;
        cam_reset = true;
        cam_x = -350;
        cam_y = -10;
        cam_z = 50;
        cam_pitch = -0.2f;
        cam_yaw = 1.1f;
        resetSliders();
    } else if (e.getSource() == preset_button_flat.get()) {
        SetupFlatPlane();
        g_reset_type = R_SQUARE;
        cam_reset = true;
        cam_x = -130;
        cam_y = -30;
        cam_z = 40;
        cam_pitch = -0.2f;
        cam_yaw = 0.785f;
        resetSliders();
    }
}


bool GuiManager::getCameraReset(float &x, float &y, float &z, float &pitch, float &yaw)
{
    if (cam_reset) {
        x = cam_x;
        y = cam_y;
        z = cam_z;
        pitch = cam_pitch;
        yaw = cam_yaw;
        cam_reset = false;
        return true;
    } else {
        return false;
    }
}

void GuiManager::resetSliders()
{
    for (int i = 0; i < int(sliders.size()); ++i) {
        Setting &s = g_settings[i];
        if (s.type == S_SLIDER || s.type == S_SLIDER_MULT_4 || s.type == S_SLIDER_INT) {
            sliders[i]->setValue(RealToSlider(s, s.value));
            SetLabelTxt(*values[i], s.value);
        } else if (s.type == S_CHECKBOX) {
            check_boxes[i]->setSelected(s.value != 0);
        }
    }
}
        
