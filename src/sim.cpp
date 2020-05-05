#include "MeiiModel.hpp"
#include <Mahi/Com.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Robo.hpp>
#include <thread>
#include <mutex>
#include <atomic>

// This is the dll that drives the model in the Unity OpenWristSim visualization
// Unity intefaces via the the C methods below, mostly just to start/stop the simulation, and read joint positions
// End-users inteface through the MelShares or at a higher level, through OpenWristSim.hpp/cpp in their own C++ code
// This was all sort of hastily put together...it could be better!

MeiiModel g_model;
std::thread g_thread;
std::mutex g_mtx;
std::atomic_bool g_stop;

using namespace mahi::com;
using namespace mahi::util;
using namespace mahi::robo;


// void stop()
// {
//     g_stop = true;
//     if (g_thread.joinable())
//         g_thread.join();
// }

// void start()
// {
//     stop();
//     g_stop = false;
//     g_thread = std::thread(simulation);
// }

// void set_torques(double tau1, double tau2, double tau3)
// {
//     std::lock_guard<std::mutex> lock(g_mtx);
//     g_model.set_torques(tau1, tau2, tau3);
// }

// void set_positions(double q1, double q2, double q3)
// {
//     std::lock_guard<std::mutex> lock(g_mtx);
//     g_model.set_positions(q1, q2, q3);
// }

// void set_velocities(double q1d, double q2d, double q3d)
// {
//     std::lock_guard<std::mutex> lock(g_mtx);
//     g_model.set_velocities(q1d, q2d, q3d);
// }

// void get_positions(double *positions)
// {
//     std::lock_guard<std::mutex> lock(g_mtx);
//     positions[0] = g_model.q1;
//     positions[1] = g_model.q2;
//     positions[2] = g_model.q3;
// }

int main(){
    MelShare ms_in("new_ms");
    std::vector<double> ms_in_data(5, 0);
    g_model.reset();
    Timer timer(hertz(500), Timer::Hybrid);
    Time t;
    Time sim_time = seconds(0);
    double kp = 500;
    double kd = 10;
    double q_ref1 = 0.1;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;

    while (true)
    {
        std::cout << "here" << std::endl;
        ms_in_data = ms_in.read_data();
        q_ref1 = ms_in_data[2];
        q_ref2 = ms_in_data[3];
        q_ref3 = ms_in_data[4]; 
        kp = ms_in_data[0];
        kd = ms_in_data[1];
        std::cout << q_ref1 << std::endl;

        double tau1 = kp * (q_ref1 - g_model.q1) - kd * g_model.q1d;
        double tau2 = kp * (q_ref2 - g_model.q2) - kd * g_model.q2d;
        double tau3 = kp * (q_ref3 - g_model.q3) - kd * g_model.q3d;
        g_model.set_torques(tau1,tau2,tau3);
        // g_model.set_torques(0.0,0.0,0.0)
        g_model.update(sim_time);
        std::cout << sim_time << ", ";
        std::cout << g_model.q1 << ", ";
        std::cout << g_model.q2 << ", ";
        std::cout << g_model.q3 << ", ";
        std::cout << g_model.q4 << ", ";
        std::cout << g_model.q5 << ", ";
        std::cout << g_model.q6 << std::endl;
        sim_time += 1_ms;
        t = timer.wait();
    }
    
    return 0;
}