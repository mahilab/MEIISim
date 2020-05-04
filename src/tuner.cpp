#define MAHI_GUI_NO_CONSOLE
#include <Mahi/Gui.hpp>
#include <Mahi/Util.hpp>
#include <thread>
#include <mutex>
#include "MeiiSim.hpp"

using namespace mahi::gui;
using namespace mahi::util;

class PdTuner : public Application {
public:
    /// Constructor
    PdTuner() : Application(500,500,"PD Tuner") 
    { 
        ImGui::DisableViewports();
        plot.x_axis.minimum = 0;   plot.x_axis.maximum = 10;
        plot.y_axis.minimum = -50; plot.y_axis.maximum = 50;
        items.resize(2);
        items[0].color = Grays::Gray50;
        items[0].data.reserve(20000);
        items[1].color = Blues::DeepSkyBlue;
        items[1].data.reserve(20000);
        control_thread = std::thread(&PdTuner::control_loop, this);
    }

    void roll_point(ImGui::PlotItem& item, Time t, double rads) {
        float tmod = fmod(t.as_seconds(), 10);
        float degs = rads * RAD2DEG;
        if (!item.data.empty() && tmod < item.data.back().x)
            item.data.clear();
        item.data.push_back({tmod, degs});
    }

    ~PdTuner() {
         stop = true;
         control_thread.join();
     }

    /// Override update from Application, called once per frame
    void update() override
    { 
        std::lock_guard<std::mutex> lock(mtx);
        ImGui::BeginFixed("PD Tuner", {0,0}, {500,500}, ImGuiWindmeiiFlags_NoTitleBar);
        if (ImGui::Checkbox("Enable", &enabled))
            meii.set_enabled(enabled);
        ImGui::SliderInt("Joint", &j, 0, 2);
        ImGui::DragDouble("Kp",&kp,0.1f,0,100);
        ImGui::DragDouble("Kd",&kd,0.01f,0,10);    
        if (ImGui::RadioButton("Square", traj.type == Waveform::Square))
            traj.type = Waveform::Square;
        ImGui::SameLine();
        if (ImGui::RadioButton("Sin", traj.type == Waveform::Sin))
            traj.type = Waveform::Sin;
        items[1].color = j == 0 ? Reds::Crimson : j == 1 ? Greens::Chartreuse : Blues::DeepSkyBlue;
        ImGui::Plot("Position", plot, items);
        ImGui::End();
    }

    void control_loop() {
        Timer timer = Timer(hertz(500)); 
        Time t = Time::Zero;    
        while(!stop) {
            {
                std::lock_guard<std::mutex> lock(mtx);
                double q  = meii.get_position(j);
                double qd = meii.get_velocity(j);
                double q_ref = traj(t);
                double tau = kp * (q_ref - q) - kd * qd;
                std::vector<double> Tau(3,0);
                Tau[j] = tau;
                meii.set_torques(Tau);
                roll_point(items[0], t, q_ref);
                roll_point(items[1], t, q);
            }
            t = timer.wait();
        }
    }

    // Member Variables
    MeiiSim meii;    
    bool enabled = false;
    int j = 0;
    double kp = 10;
    double kd = 1;
    Waveform traj = Waveform(Waveform::Square, seconds(2), 30 * DEG2RAD);
    ImGui::PlotInterface plot;
    std::vector<ImGui::PlotItem> items;
    std::thread control_thread;
    std::atomic_bool stop = false;
    std::mutex mtx;
};

int main(int argc, char const *argv[])
{
    PdTuner tuner;
    tuner.run();
    return 0;
}