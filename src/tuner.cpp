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

struct ScrollingData {
    int MaxSize = 1000;
    int Offset  = 0;
    ImVector<ImVec2> Data;
    ScrollingData() { Data.reserve(MaxSize); }
    void AddPoint(float x, float y) {
        if (Data.size() < MaxSize)
            Data.push_back(ImVec2(x,y));
        else {
            Data[Offset] = ImVec2(x,y);
            Offset =  (Offset + 1) % MaxSize;
        }
    }
};

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

    ~PdTuner() {
         stop = true;
     }

    /// Override update from Application, called once per frame
    void update() override
    {
        std::vector<double> ms_in_data = ms_in.read_data();
        ImGui::Begin("PD Tuner");
        ImGui::DragDouble("Kp",&kp,0,0,1000);
        ImGui::DragDouble("Kd",&kd,0,0,100);
        ImGui::DragDouble("Khard",&k_hard,0,0,20000);
        ImGui::DragDouble("Bhard",&b_hard,0,0,2000);
        ImGui::DragDouble("q_ref 1",&q_ref1,0.0001f,0.03,0.15,"%.4f");
        ImGui::DragDouble("q_ref 2",&q_ref2,0.0001f,0.03,0.15,"%.4f");    
        ImGui::DragDouble("q_ref 3",&q_ref3,0.0001f,0.03,0.15,"%.4f");    
        ImGui::DragDouble("q1",&q1,0.0001f,0.03,0.15,"%.4f");
        ImGui::DragDouble("q2",&q2,0.0001f,0.03,0.15,"%.4f");    
        ImGui::DragDouble("q3",&q3,0.0001f,0.03,0.15,"%.4f");    

        if(!ms_in_data.empty()) {
            sim_rate.AddPoint(t,ms_in_data[0]);
            tau1.AddPoint(t,ms_in_data[1]);
            tau2.AddPoint(t,ms_in_data[2]);
            tau3.AddPoint(t,ms_in_data[3]);
            q1 = ms_in_data[4];
            q2 = ms_in_data[5];
            q3 = ms_in_data[6];
        }
        ImGui::SetNextPlotRangeX(t - 10, t, ImGuiCond_Always);
        ImGui::SetNextPlotRangeY(1500,3500);
        ImGui::BeginPlot("##Sim Time", "Time (s)", "Sim Rate (us)", {-1,200}, ImPlotFlags_Default, rt_axis, rt_axis);
        ImGui::Plot("Sim Time", &sim_rate.Data[0].x, &sim_rate.Data[0].y, sim_rate.Data.size(), sim_rate.Offset, 2 * sizeof(float));
        ImGui::EndPlot();

        ImGui::SetNextPlotRangeX(t - 10, t, ImGuiCond_Always);
        ImGui::SetNextPlotRangeY(-200,200);
        ImGui::BeginPlot("##Forces", "Time (s)", "Force (N)", {-1,400}, ImPlotFlags_Default, rt_axis, rt_axis);
        ImGui::Plot("Tau 1", &tau1.Data[0].x, &tau1.Data[0].y, tau1.Data.size(), tau1.Offset, 2 * sizeof(float));
        ImGui::Plot("Tau 2", &tau2.Data[0].x, &tau2.Data[0].y, tau2.Data.size(), tau2.Offset, 2 * sizeof(float));
        ImGui::Plot("Tau 3", &tau3.Data[0].x, &tau3.Data[0].y, tau3.Data.size(), tau3.Offset, 2 * sizeof(float));
        ImGui::EndPlot();

        t += ImGui::GetIO().DeltaTime;
        ImGui::End();
        std::vector<double> data = {kp, kd, q_ref1, q_ref2, q_ref3, k_hard, b_hard};
        ms_out.write_data(data);
    }

    // Member Variables
    int rt_axis = ImAxisFlags_Default & ~ImAxisFlags_TickLabels;
    ScrollingData sim_rate;
    ScrollingData tau1;
    ScrollingData tau2;
    ScrollingData tau3;
    float t = 0;
    bool enabled = false;
    double k_hard = 100;
    double b_hard = 10;
    double kp = 400;
    double kd = 8;
    double q_ref1 = 0.1;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    double q1 = 0.1;
    double q2 = 0.1;
    double q3 = 0.1;
    MelShare ms_out = MelShare("sim1");
    MelShare ms_in = MelShare("sim2");
    std::atomic_bool stop = false;
    std::mutex mtx;
};

int main(int argc, char const *argv[])
{
    PdTuner tuner;
    tuner.run();
    return 0;
}