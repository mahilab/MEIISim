#pragma once
#include <Mahi/Util/Timing/Time.hpp>
#include <Mahi/Util/Math/Constants.hpp>
#include <Mahi/Util/Math/Integrator.hpp>
#include <Mahi/Robo/Control/Limiter.hpp>

/// Dynamic Model of the OpenWrist
class MeiiModel {
public:

    MeiiModel();

    void update(mahi::util::Time t);
    void set_torques(double tau1, double tau2, double tau3);
    void set_positions(double _q1, double _q2, double _q3, double _q4, double _q5, double _q6);
    void set_velocities(double _q1d, double _q2d, double _q3d, double _q4d, double _q5d, double _q6d);
    void reset();
    void calc_dependent_joint_values();

public:

    // Joint torques [Nm]
    double tau1, tau2, tau3;
    // Joint Positions [m]
    double q1, q2, q3;
    // Joint Velocities [m/s]
    double q1d, q2d, q3d;
    // Joint Accelerations [m/s^2]
    double q1dd, q2dd, q3dd;

    // Dependent joint positions [rad]
    // q4 -> theta1
    // q5 -> theta2
    // q6 -> theta3
    // q7 -> Px
    // q8 -> Py
    // q9 -> Pz
    // q10 -> alpha (Rx)
    // q11 -> beta (Ry)
    // q12 -> gamma (Rz)
    double q4, q5, q6, q7, q8, q9, q10, q11, q12;
    // Dependent joint positions [rad/s]
    double q4d, q5d, q6d;

    // Hardstops
    const double q1min = 0.050;
    const double q2min = 0.050;
    const double q3min = 0.050;
    const double q1max = 0.133;
    const double q2max = 0.133;
    const double q3max = 0.133;

    const double Khard = 20000; // hardstop stiffness
    const double Bhard = 1000;  // hardstop damping

    // // Joint Mass [kg]
    // const double m1 = 1.79265300000000;
    // const double m2 = 0.891430000000000;
    // const double m3 = 0.192063000000000;

    // // Joint Moments of Inertia [kg*m^2]
    // const double Ic1xx = 0.00530900000000000;  
    // const double Ic1xy = -0.000758000000000000; 
    // const double Ic1xz = -0.000687000000000000; 
    // const double Ic1yy = 0.00899500000000000; 
    // const double Ic1yz = -0.000368000000000000; 
    // const double Ic1zz = 0.0103780000000000; 
    // const double Ic2xx = 0.00375400000000000; 
    // const double Ic2xy = -0.000527000000000000; 
    // const double Ic2xz = 7.10000000000000e-05; 
    // const double Ic2yy = 0.00151900000000000; 
    // const double Ic2yz = -0.000324000000000000; 
    // const double Ic2zz = 0.00451600000000000; 
    // const double Ic3xx = 0.000419000000000000; 
    // const double Ic3xy = -2.60000000000000e-05; 
    // const double Ic3xz = -0.000140000000000000; 
    // const double Ic3yy = 0.000470000000000000; 
    // const double Ic3yz = -2.90000000000000e-05; 
    // const double Ic3zz = 0.000287000000000000; 

    // // Joint Center of Mass [m]
    // const double Pc1x = 0.0224710000000000; 
    // const double Pc1y = 0.0404220000000000;
    // const double Pc1z = 0.115783000000000;
    // const double Pc2x = -0.00945400000000000;
    // const double Pc2y = -0.0271490000000000;
    // const double Pc2z = -0.0781620000000000;
    // const double Pc3x = 0.0493680000000000;
    // const double Pc3y = -0.0211170000000000;
    // const double Pc3z = 0.0607280000000000;

    // // Motor Rotor Inertia [kg-m^2]
    // const double Jm1 = 1.37000000000000e-05;
    // const double Jm2 = 1.37000000000000e-05;
    // const double Jm3 = 3.47000000000000e-06;

    const double R = 0.1044956; // m
    const double r = 0.05288174521; // m
    const double a5 = 0.0268986; // m
    const double a6 = 0.0272820; // m
    const double a56 = -(a5-a6); // m
    const double alpha_5 = 0.094516665; // rad
    const double alpha_13 = 5*mahi::util::PI/180; // rad

    /// Motor continous and max torque limits [Nm]
    const double tau1_mot_cont = 0.187; 
    const double tau1_mot_max  = 2.560;
    const double tau2_mot_cont = 0.187; 
    const double tau2_mot_max  = 2.560;
    const double tau3_mot_cont = 0.0897; 
    const double tau3_mot_max  = 1.050;

    // Transmission Rations [rad/m]
    const double eta1 = 0.23*0.0254;
    const double eta2 = 0.23*0.0254;
    const double eta3 = 0.23*0.0254;

    // Damping Coefficients [Nm*s/rad]
    const double b1 = 0.5 * 0.0252;  
    const double b2 = 0.5 * 0.0019;  
    const double b3 = 0.5 * 0.0029; 

    // Kinetic Friction [Nm]
    const double fk1 = 0.5 * 0.1891; 
    const double fk2 = 0.5 * 0.0541;
    const double fk3 = 0.5 * 0.1339; 

    // Gravity Constant [m/s^2]
    const double g = 9.80665;

private:
    // torque limiters
    mahi::robo::Limiter lim1, lim2, lim3;

    // Integrators
    mahi::util::Integrator q1dd_q1d, q2dd_q2d, q3dd_q3d, q1d_q1, q2d_q2, q3d_q3;
};