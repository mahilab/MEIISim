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
    MelShare ms_gains("gains");
    MelShare ms_refs("refs");
    MelShare ms_times_out("times");
    MelShare ms_qs_out("qs");
    double kp = 400;
    double kd = 15;
    double kpf = 20;
    double kdf = 2;
    double kpe = 30;
    double kde = 3;
    double q_ref1 = 0.0;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    double q_ref4 = 0.1;
    double q_ref5 = 0.5;
    double q1, q2, q3, q4, q5;
    double tau1, tau2, tau3, tau4, tau5;
    double k_hard1 = 200;
    double b_hard1 = 10;
    double threadpooling = true;
    Timer timer(hertz(1000), Timer::Hybrid);
    Time t;
    Time t_last;
    Time sim_time = 0_ms;
    double calc_time = 0;
    double setup_time = 0;
    double comp_time = 0;
    int n_threads = 12;
    int n_threads_last = n_threads;
    while (!g_stop)
    {
        auto ms_gain_data = ms_gains.read_data();
        auto ms_ref_data = ms_refs.read_data();
        if (!ms_gain_data.empty()){
            kp = ms_gain_data[0];
            kd = ms_gain_data[1];
            kpe = ms_gain_data[2];
            kde = ms_gain_data[3];
            kpf = ms_gain_data[4];
            kdf = ms_gain_data[5];
            k_hard1 = ms_gain_data[6];
            b_hard1 = ms_gain_data[7];

            q_ref1 = ms_ref_data[0];
            q_ref2 = ms_ref_data[1];
            q_ref3 = ms_ref_data[2]; 
            q_ref4 = ms_ref_data[3]; 
            q_ref5 = ms_ref_data[4]; 
            
            threadpooling = ms_ref_data[5];
            n_threads = int(ms_ref_data[6]);
        }
        {
            std::lock_guard<std::mutex> lock(g_mtx);
            g_model.threadpool = (threadpooling > 0.5) ? true : false;
            if (n_threads != n_threads_last) g_model.p.resize(n_threads);
            n_threads_last = n_threads;
            tau1 = kpe * (q_ref1 - g_model.q1) - kde * g_model.q1d;
            tau2 = kpf * (q_ref2 - g_model.q2) - kdf * g_model.q2d;
            tau3 = kp  * (q_ref3 - g_model.q3) - kd  * g_model.q3d;
            tau4 = kp  * (q_ref4 - g_model.q4) - kd  * g_model.q4d;
            tau5 = kp  * (q_ref5 - g_model.q5) - kd  * g_model.q5d;
            g_model.Khard1 = k_hard1;
            g_model.Bhard1 = b_hard1;
            g_model.set_torques(tau1,tau2,tau3,tau4,tau5);
            g_model.update(sim_time);
            q1 = g_model.q1;
            q2 = g_model.q2;
            q3 = g_model.q3;
            q4 = g_model.q4;
            q5 = g_model.q5;
            calc_time = g_model.mat_calc_time;
            setup_time = g_model.setup_time;
            comp_time = g_model.comp_time;
        }
        sim_time += 1_ms;
        ms_times_out.write_data({double((t-t_last).as_microseconds()),calc_time,setup_time,comp_time});
        ms_qs_out.write_data({tau1,tau2,tau3,tau4,tau5,q1,q2,q3,q4,q5});
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

EXPORT void set_torques(double tau1, double tau2, double tau3, double tau4, double tau5)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_torques(tau1, tau2, tau3, tau4, tau5);
}

EXPORT void set_positions(double q1, double q2, double q3, double q4, double q5, double q6, double q7, double q8)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_positions(q1, q2, q3, q4, q5, q6, q7, q8);
}

EXPORT void set_velocities(double q1d, double q2d, double q3d, double q4d, double q5d, double q6d, double q7d, double q8d)
{
    std::lock_guard<std::mutex> lock(g_mtx);
    g_model.set_velocities(q1d, q2d, q3d, q4d, q5d, q6d, q7d, q8d);
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
    positions[12] = g_model.q13;
    positions[13] = g_model.q14;
}