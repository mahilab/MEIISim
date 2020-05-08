#define MAHI_GUI_NO_CONSOLE
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

class SimTuner : public Application {
public:
    /// Constructor
    SimTuner() : Application(500,500,"Sim Tuner"),
    calc_times_1s(60)
    { 
        sim_rate.AddPoint(t,0);
        des_rate.AddPoint(t,0);
        tau1.AddPoint(t,0);
        tau2.AddPoint(t,0);
        tau3.AddPoint(t,0);
        calcTime.AddPoint(t,0);
        setupTime.AddPoint(t,0);
        compTime.AddPoint(t,0);
    }

    ~SimTuner() {
         stop = true;
     }

    /// Override update from Application, called once per frame
    void update() override
    {
        std::vector<double> ms_in_data = ms_in.read_data();
        if(!ms_in_data.empty()) {
            sim_rate.AddPoint(t,ms_in_data[0]);
            des_rate.AddPoint(t,1000.0f);
            tau1.AddPoint(t,ms_in_data[1]);
            tau2.AddPoint(t,ms_in_data[2]);
            tau3.AddPoint(t,ms_in_data[3]);
            q1 = ms_in_data[4];
            q2 = ms_in_data[5];
            q3 = ms_in_data[6];
            calcTime.AddPoint(t,ms_in_data[7]);
            calc_times_1s.push_back(ms_in_data[7]);
            setupTime.AddPoint(t,ms_in_data[8]);
            compTime.AddPoint(t,ms_in_data[9]);
        }
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
        ImGui::Checkbox("Threadpooled", &tp);    
        ImGui::SliderInt("Number of Threads", &nthread, 1, 20);
        double mean_calc_time = mean(calc_times_1s.get_vector());
        ImGui::InputDouble("Mean Calc Time",&mean_calc_time);

        ImGui::SetNextPlotRangeX(t - 10, t, ImGuiCond_Always);
        ImGui::SetNextPlotRangeY(0,1500);
        if(ImGui::BeginPlot("##Sim Time", "Time (s)", "Sim Rate (us)", {-1,300}, ImPlotFlags_Default, rt_axis, rt_axis)){
            ImGui::Plot("Sim Loop Time", &sim_rate.Data[0].x, &sim_rate.Data[0].y, sim_rate.Data.size(), sim_rate.Offset, 2 * sizeof(float));
            ImGui::Plot("Desired Sim Loop Time", &des_rate.Data[0].x, &des_rate.Data[0].y, des_rate.Data.size(), des_rate.Offset, 2 * sizeof(float));
            ImGui::Plot("Total Calc Time", &calcTime.Data[0].x, &calcTime.Data[0].y, calcTime.Data.size(), calcTime.Offset, 2 * sizeof(float));
            ImGui::Plot("Setup Time (thread only)", &setupTime.Data[0].x, &setupTime.Data[0].y, setupTime.Data.size(), setupTime.Offset, 2 * sizeof(float));
            ImGui::Plot("Comp Time (thread only)", &compTime.Data[0].x, &compTime.Data[0].y, compTime.Data.size(), compTime.Offset, 2 * sizeof(float));
            ImGui::EndPlot();
        }

        ImGui::SetNextPlotRangeX(t - 10, t, ImGuiCond_Always);
        ImGui::SetNextPlotRangeY(-10,10);
        if(ImGui::BeginPlot("##Forces", "Time (s)", "Force (N)", {-1,300}, ImPlotFlags_Default, rt_axis, rt_axis)){
            ImGui::Plot("Force 1", &tau1.Data[0].x, &tau1.Data[0].y, tau1.Data.size(), tau1.Offset, 2 * sizeof(float));
            ImGui::Plot("Force 2", &tau2.Data[0].x, &tau2.Data[0].y, tau2.Data.size(), tau2.Offset, 2 * sizeof(float));
            ImGui::Plot("Force 3", &tau3.Data[0].x, &tau3.Data[0].y, tau3.Data.size(), tau3.Offset, 2 * sizeof(float));
            ImGui::EndPlot();
        }

        t += ImGui::GetIO().DeltaTime;
        ImGui::End();
        double tp_double = tp ? 1.0 : 0.0;
        std::vector<double> data = {kp, kd, q_ref1, q_ref2, q_ref3, k_hard, b_hard, tp_double, double(nthread)};
        ms_out.write_data(data);
    }

    // Member Variables
    int rt_axis = ImAxisFlags_Default & ~ImAxisFlags_TickLabels;
    ScrollingData sim_rate;
    ScrollingData des_rate;
    ScrollingData tau1;
    ScrollingData tau2;
    ScrollingData tau3;
    ScrollingData calcTime;
    ScrollingData setupTime;
    ScrollingData compTime;
    float t = 0;
    bool tp = false;
    bool enabled = false;
    double k_hard = 20000;
    double b_hard = 100;
    double kp = 600;
    double kd = 15;
    double q_ref1 = 0.1;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    int nthread = 5;
    double q1 = 0.1;
    double q2 = 0.1;
    double q3 = 0.1;
    MelShare ms_out = MelShare("sim1");
    MelShare ms_in = MelShare("sim2");
    std::atomic_bool stop = false;
    std::mutex mtx;
    RingBuffer<double> calc_times_1s;
};

int main(int argc, char const *argv[])
{
    SimTuner tuner;
    tuner.run();
    return 0;
}