#pragma once
#include <Mahi/Util/Timing/Time.hpp>
#include <Mahi/Util/Math/Constants.hpp>
#include <Mahi/Util/Math/Integrator.hpp>
#include <Mahi/Robo/Control/Limiter.hpp>
#include <Eigen/Dense>
#include <ctpl_stl.h>

/// Dynamic Model of the OpenWrist
class MeiiModel {
public:

    MeiiModel();

    void update(mahi::util::Time t);
    void set_torques(double tau1, double tau2, double tau3, double tau4, double tau5);
    void set_positions(double _q1, double _q2, double _q3, double _q4, double _q5, double _q6, double _q7, double _q8);
    void set_velocities(double _q1d, double _q2d, double _q3d, double _q4d, double _q5d, double _q6d, double _q7d, double _q8d);
    void reset();
    void calc_dependent_joint_values();

public:

    bool threadpool = true;

    // Joint torques [Nm]
    double tau1, tau2, tau3, tau4, tau5;
    // Joint Positions [m]
    double   q1,   q2,   q3,   q4,   q5;
    // Joint Velocities [m/s]
    double  q1d,  q2d,  q3d,  q4d,  q5d;
    // Joint Accelerations [m/s^2]
    double q1dd, q2dd, q3dd, q4dd, q5dd;

    // Dependent joint positions [rad]
    // q5 -> theta1
    // q6 -> theta2
    // q7 -> theta3
    // q8 -> Px
    // q9 -> Py
    // q10 -> Pz
    // q11 -> alpha (Rx)
    // q12 -> beta (Ry)
    // q13 -> gamma (Rz)
    double q6, q7, q8, q9, q10, q11, q12, q13, q14;
    // Dependent joint velocities [rad/s]
    double q6d, q7d, q8d;

    // Hardstops
    const double q1min = -91.5 * mahi::util::DEG2RAD;
    const double q1max = 3.0 * mahi::util::DEG2RAD;
    const double q2min = -99 * mahi::util::DEG2RAD;
    const double q2max = 108 * mahi::util::DEG2RAD;
    const double q3min = 0.050;
    const double q4min = 0.050;
    const double q5min = 0.050;
    const double q3max = 0.1305;
    const double q4max = 0.1305;
    const double q5max = 0.1305;

    //// Hardstops if you want to see the effect of gravity without instablities
    // const double q3min = 0.09;
    // const double q4min = 0.09;
    // const double q5min = 0.09;
    // const double q3max = 0.110;
    // const double q4max = 0.110;
    // const double q5max = 0.110;

    double Khard = 20000; // hardstop stiffness
    double Bhard = 100;  // hardstop damping
    double Khard1 = 2000; // hardstop stiffness
    double Bhard1 = 100;  // hardstop damping
    double mat_calc_time = 0;
    double setup_time = 0;
    double comp_time = 0;

    const double R = 0.1044956; // m
    const double r = 0.05288174521; // m
    const double a4 = 0.159385;
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
    const double tau4_mot_cont = 0.187; 
    const double tau4_mot_max  = 2.560;
    const double tau5_mot_cont = 0.0897; 
    const double tau5_mot_max  = 1.050;

    // Transmission Rations [rad/m]
    const double eta1 = 0.23*0.0254;
    const double eta2 = 0.23*0.0254;
    const double eta3 = 0.23*0.0254;
    const double eta4 = 0.23*0.0254;
    const double eta5 = 0.23*0.0254;

    // Damping Coefficients [Nm*s/rad]
    const double b1 = 0.5 * 0.0252;  
    const double b2 = 0.5 * 0.0019;  
    const double b3 = 0.5 * 0.0029; 
    const double b4 = 0.5 * 0.0019;  
    const double b5 = 0.5 * 0.0029; 

    // Kinetic Friction [Nm]
    const double fk1 = 0.5 * 0.1891;
    const double fk2 = 0.5 * 0.0541;
    const double fk3 = 0.5 * 0.1339;
    const double fk4 = 0.5 * 0.1339;
    const double fk5 = 0.5 * 0.1339;

    // Gravity Constant [m/s^2]
    const double g = 9.80665;

    ctpl::thread_pool p;

private:
    // torque limiters
    mahi::robo::Limiter lim1, lim2, lim3, lim4, lim5;

    // Integrators
    mahi::util::Integrator q1dd_q1d, q2dd_q2d, q3dd_q3d, q4dd_q4d, q5dd_q5d, q1d_q1, q2d_q2, q3d_q3, q4d_q4, q5d_q5;

    std::future<double> v00, v01, v02, v03, v04, v05, v06, v07, v10, v11, v12, v13, v14, v15, v16, v17, v20, v21, v22, v23, v24, v25, v26, v27, v30, v31, v32, v33, v34, v35, v36, v37, v40, v41, v42, v43, v44, v45, v46, v47, v50, v51, v52, v53, v54, v55, v56, v57, v60, v61, v62, v63, v64, v65, v66, v67, v70, v71, v72, v73, v74, v75, v76, v77;

    std::future<double> m00, m01, m02, m03, m04, m05, m06, m07, m10, m11, m12, m13, m14, m15, m16, m17, m20, m21, m22, m23, m24, m25, m26, m27, m30, m31, m32, m33, m34, m35, m36, m37, m40, m41, m42, m43, m44, m45, m46, m47, m50, m51, m52, m53, m54, m55, m56, m57, m60, m61, m62, m63, m64, m65, m66, m67, m70, m71, m72, m73, m74, m75, m76, m77;

    std::future<double> g0, g1, g2, g3, g4, g5, g6, g7;

    Eigen::MatrixXd V;
    Eigen::MatrixXd M;
    Eigen::VectorXd G;
    Eigen::VectorXd Tau;
    Eigen::VectorXd B;
    Eigen::VectorXd Fk;

    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd x;

};