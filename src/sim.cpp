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
    g_model.reset();
    Timer timer(hertz(500), Timer::Hybrid);
    Time t;
    double q_ref = 0.1;
    double kp = 5000;
    double kd = 10;
    while (true)
    {
        
        double tau1 = kp * (q_ref - g_model.q1) - kd * g_model.q1d;
        double tau2 = kp * (q_ref - g_model.q2) - kd * g_model.q2d;
        double tau3 = kp * (q_ref - g_model.q3) - kd * g_model.q3d;
        g_model.set_torques(tau1,tau2,tau3);
        // g_model.set_torques(0.0,0.0,0.0);
        g_model.update(t);
        std::cout << t << ", ";
        std::cout << g_model.q1 << ", ";
        std::cout << g_model.q2 << ", ";
        std::cout << g_model.q3 << ", ";
        std::cout << g_model.q4 << ", ";
        std::cout << g_model.q5 << ", ";
        std::cout << g_model.q6 << std::endl;
        t = timer.wait();
    }
    
    return 0;
}