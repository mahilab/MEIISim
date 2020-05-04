#pragma once
#include <Mahi/Com/MelShare.hpp>
#include <Mahi/Robo/Control/PdController.hpp>
#include <vector>

/// High-level interface to interact with the Unity MeiiSim visualization
class MeiiSim {
public:
    MeiiSim();
    ~MeiiSim();
    void set_enabled(bool enabled);
    double get_position(int i);
    double get_velocity(int i);
    double get_acceleration(int i);
    double get_torque(int i);
    void set_torque(int i, double tau);
    std::vector<double> get_positions();
    std::vector<double> get_velocities();
    std::vector<double> get_accelerations();
    std::vector<double> get_torques();
    void set_torques(const std::vector<double>& Tau);
    void calibrate();
public:
    std::vector<mahi::robo::PdController> pd = std::vector<mahi::robo::PdController>{
        mahi::robo::PdController(25, 1.15), // joint 0 ( Nm/rad , Nm-s/rad )
        mahi::robo::PdController(20, 1.00), // joint 1 ( Nm/rad , Nm-s/rad )
        mahi::robo::PdController(20, 0.25)  // joint 2 ( Nm/rad , Nm-s/rad )
    };
private:
    bool enabled = false;
    mahi::com::MelShare ms_in, ms_out;
    std::vector<double> ms_in_data, ms_out_data;
};