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

#define EXPORT extern "C" __declspec(dllexport)

MeiiModel g_model;
std::thread g_thread;
std::mutex g_mtx;
std::atomic_bool g_stop;

using namespace mahi::com;
using namespace mahi::util;
using namespace mahi::robo;

void simulation()
{
    MelShare ms_in("meii_sim_kp_kd_ref");
    MelShare ms_out("meii_sim_rate_tau");
    std::vector<double> ms_in_data(5, 0);
    double kp = 500;
    double kd = 10;
    double q_ref1 = 0.1;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    double tau1 = 0;
    double tau2 = 0;
    double tau3 = 0;
    Timer timer(hertz(500), Timer::Hybrid);
    Time t;
    Time t_last;
    Time sim_time = 0_ms;
    while (!g_stop)
    {
        ms_in_data = ms_in.read_data();
        if (!ms_in_data.empty()){
            q_ref1 = ms_in_data[2];
            q_ref2 = ms_in_data[3];
            q_ref3 = ms_in_data[4]; 
            kp = ms_in_data[0];
            kd = ms_in_data[1];
        }
        {
            std::lock_guard<std::mutex> lock(g_mtx);
            tau1 = kp * (q_ref1 - g_model.q1) - kd * g_model.q1d;
            tau2 = kp * (q_ref2 - g_model.q2) - kd * g_model.q2d;
            tau3 = kp * (q_ref3 - g_model.q3) - kd * g_model.q3d;
            g_model.set_torques(tau1,tau2,tau3);
            g_model.update(sim_time);
        }
        sim_time += 1_ms;
        ms_out.write_data({double((t-t_last).as_microseconds()),tau1,tau2,tau3});
        t_last = t;
        t = timer.wait();
    }
};

EXPORT void stop()
{
    g_stop = true;
    if (g_thread.joinable())
        g_thread.join();
}

EXPORT void start()
{
    stop();
    g_model.reset();
    g_stop = false;
    g_thread = std::thread(simulation);
}

EXPORT void set_torques(double tau1, double tau2, double tau3)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_torques(tau1, tau2, tau3);
}

EXPORT void set_positions(double q1, double q2, double q3, double q4, double q5, double q6)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_positions(q1, q2, q3, q4, q5, q6);
}

EXPORT void set_velocities(double q1d, double q2d, double q3d, double q4d, double q5d, double q6d)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_velocities(q1d, q2d, q3d, q4d, q5d, q6d);
}

EXPORT void get_positions(double *positions)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    positions[0] = g_model.q1;
    positions[1] = g_model.q2;
    positions[2] = g_model.q3;
    positions[3] = g_model.q4;
    positions[4] = g_model.q5;
    positions[5] = g_model.q6;
    positions[6] = g_model.q7;
    positions[7] = g_model.q8;
    positions[8] = g_model.q9;
    positions[9] = g_model.q10;
    positions[10] = g_model.q11;
    positions[11] = g_model.q12;
}