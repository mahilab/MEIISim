#include "MeiiModel.hpp"
#include <Mahi/Com.hpp>
#include <Mahi/Util.hpp>
#include <Mahi/Robo.hpp>
#include <thread>
#include <mutex>
#include <atomic>

MeiiModel g_model;
std::atomic_bool g_stop;

using namespace mahi::com;
using namespace mahi::util;
using namespace mahi::robo;

int main(){
    double kp = 600;
    double kd = 15;
    double kpf = 20;
    double kdf = 0.5;
    double q_ref1 = 0.0;
    double q_ref2 = 0.1;
    double q_ref3 = 0.1;
    double q_ref4 = 0.1;
    double q1, q2, q3, q4;
    double tau1 = 0;
    double tau2 = 0;
    double tau3 = 0;
    double tau4 = 0;
    Time sim_time = 0_s;
    g_model.threadpool = true;
    while (!g_stop) {
        tau1 = 0;
        tau2 = kp * (q_ref2 - g_model.q2) - kd * g_model.q2d;
        tau3 = kp * (q_ref3 - g_model.q3) - kd * g_model.q3d;
        tau4 = kp * (q_ref4 - g_model.q4) - kd * g_model.q4d;
        g_model.set_torques(tau1, tau2, tau3, tau4);
        g_model.update(sim_time);
        std::cout << sim_time.as_milliseconds() << std::endl;
        std::cout << g_model.q1  << ", ";
        std::cout << g_model.q2  << ", ";
        std::cout << g_model.q3  << ", ";
        std::cout << g_model.q4  << ", ";
        std::cout << g_model.q5  << ", ";
        std::cout << g_model.q6  << ", ";
        std::cout << g_model.q7  << ", ";
        std::cout << g_model.q8  << ", ";
        std::cout << g_model.q9  << ", ";
        std::cout << g_model.q10 << ", ";
        std::cout << g_model.q11 << ", ";
        std::cout << g_model.q12 << ", ";
        std::cout << g_model.q13 << std::endl;
        sleep(1_ms);
        sim_time += 1_ms;
        // std::cin.get();
    }
    
    return 0;
}