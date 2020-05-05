// #define MAHI_GUI_NO_CONSOLE
#include <Mahi/Gui.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Com.hpp>
#include <thread>
#include <mutex>
#include <atomic>

using namespace mahi::gui;
using namespace mahi::util;
using namespace mahi::com;

class PdTuner : public Application {
public:
    /// Constructor
    PdTuner() : Application(500,500,"PD Tuner") 
    { 
        // ImGui::DisableViewports();
        // plot.x_axis.minimum = 0;   plot.x_axis.maximum = 10;
        // plot.y_axis.minimum = -50; plot.y_axis.maximum = 50;
        // items.resize(2);
        // items[0].color = Grays::Gray50;
        // items[0].data.reserve(20000);
        // items[1].color = Blues::DeepSkyBlue;
        // items[1].data.reserve(20000);
        // control_thread = std::thread(&PdTuner::control_loop, this);
    }

    // void roll_point(ImGui::PlotItem& item, Time t, double rads) {
    //     float tmod = fmod(t.as_seconds(), 10);
    //     float degs = rads * RAD2DEG;
    //     if (!item.data.empty() && tmod < item.data.back().x)
    //         item.data.clear();
    //     item.data.push_back({tmod, degs});
    // }

    ~PdTuner() {
         stop = true;
     }

    /// Override update from Application, called once per frame
    void update() override
    { 
        ImGui::Begin("PD Tuner");
        // if (ImGui::Checkbox("Enable", &enabled))
        //     meii.set_enabled(enabled);
        // ImGui::SliderInt("Joint", &j, 0, 2);
        ImGui::DragDouble("Kp",&kp,0,0,1000);
        ImGui::DragDouble("Kd",&kd,0,0,100);
        ImGui::DragDouble("q_ref 1",&q_ref1,0.0001f,0.03,0.15);
        ImGui::DragDouble("q_ref 2",&q_ref2,0.0001f,0.03,0.15);    
        ImGui::DragDouble("q_ref 3",&q_ref3,0.0001f,0.03,0.15);    
        // if (ImGui::RadioButton("Square", traj.type == Waveform::Square))
            // traj.type = Waveform::Square;
        // ImGui::SameLine();
        // if (ImGui::RadioButton("Sin", traj.type == Waveform::Sin))
            // traj.type = Waveform::Sin;
        // items[1].color = j == 0 ? Reds::Crimson : j == 1 ? Greens::Chartreuse : Blues::DeepSkyBlue;
        // ImGui::Plot("Position", plot, items);
        std::vector<double> data = {kp, kd, q_ref1, q_ref2, q_ref3};
        ms_out.write_data(data);
        // std::cout << kp << ", " << kd << ", " << q_ref1 << ", " << q_ref2 << ", " << q_ref3 << std::endl;
        ImGui::End();
    }

    // Member Variables
    bool enabled = false;
    double kp = 400;
    double kd = 8;
    double q_ref1 = 0.1;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    MelShare ms_out = MelShare("meii_sim");
    std::atomic_bool stop = false;
    std::mutex mtx;
};

int main(int argc, char const *argv[])
{
    PdTuner tuner;
    tuner.run();
    return 0;
}