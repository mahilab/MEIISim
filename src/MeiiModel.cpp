#include "MeiiModel.hpp"
#include <Eigen/Dense>

using Eigen::Matrix3d;
using Eigen::Vector3d;

using namespace mahi::util;

MeiiModel::OpenWristModel() : 
{
    reset();
}

void MeiiModel::update(Time t)
{
    // limit torques
    tau1 = lim1.limit(tau1);
    tau2 = lim2.limit(tau2);
    tau3 = lim3.limit(tau3);

    // get list of qs and q_dots
    std::vector<double> qs = {q1 q2 q3 q4 q5 q6 q1d q2d q3d q4d q5d q6d}

    // Mass matrix
    Matrix3d M;
    M(0,0) = get_M11(qs);
    M(0,1) = get_M12(qs);
    M(0,2) = get_M13(qs);
    M(1,0) = get_M21(qs);
    M(1,1) = get_M22(qs);
    M(1,2) = get_M23(qs);
    M(2,0) = get_M31(qs);
    M(2,1) = get_M32(qs);
    M(2,2) = get_M33(qs);

    // Reflected motor rotor inertia
    // Matrix3d M_mot =  Matrix3d::Zero();
    // M_mot(0,0) = Jm1*eta1*eta1;
    // M_mot(1,1) = Jm2*eta2*eta2;
    // M_mot(2,2) = Jm3*eta3*eta3;

    // Coriolis vector
    Matrix3d V;
    V(0,0) = get_V11(qs);
    V(0,1) = get_V12(qs);
    V(0,2) = get_V13(qs);
    V(1,0) = get_V21(qs);
    V(1,1) = get_V22(qs);
    V(1,2) = get_V23(qs);
    V(2,0) = get_V31(qs);
    V(2,1) = get_V32(qs);
    V(2,2) = get_V33(qs);
   
    // Gravity vector
    Vector3d G;
    G[0] = get_G1(qs);
    G[1] = get_G1(qs);
    G[2] = get_G1(qs);

    // Damping
    // Vector3d B;
    // B[0] = b1*q1d;
    // B[1] = b2*q2d;
    // B[2] = b3*q3d;

    // Kinetic friction
    // Vector3d Fk;
    // Fk[0] = fk1*std::tanh(q1d*10);
    // Fk[1] = fk2*std::tanh(q2d*10);
    // Fk[2] = fk3*std::tanh(q3d*10);    

    // Torque vector
    Vector3d Tau;
    Tau[0] = tau1;// + hardstop_torque(q1,q1d,q1min,q1max,Khard,Bhard);
    Tau[1] = tau2;// + hardstop_torque(q2,q2d,q2min,q2max,Khard,Bhard);
    Tau[2] = tau3;// + hardstop_torque(q3,q3d,q3min,q3max,Khard,Bhard);

    // Solved for accelerations
    // 1) Tau = (M + M_mot) * Qdd + V + G + B + Fk
    // 2) (M + M_mot) * Qdd = Tau - V -G - B - Fk
    // 3) A             x   = b
    Matrix3d A = M;// + M_mot;
    Vector3d b = Tau - V - G;// - B - Fk;
    Vector3d x = A.householderQr().solve(b);

    q1dd = x[0];
    q2dd = x[1];
    q3dd = x[2];

    // integrate acclerations to find velocities
    q1d = q1dd_q1d.update(q1dd, t);
    q2d = q2dd_q2d.update(q2dd, t);
    q3d = q3dd_q3d.update(q3dd, t);

    // integrate velocities to find positions
    q1 = q1d_q1.update(q1d, t);
    q2 = q2d_q2.update(q2d, t);
    q3 = q3d_q3.update(q3d, t);

}

void OpenWristModel::set_torques(double _tau1, double _tau2, double _tau3) {
    tau1 = _tau1;
    tau2 = _tau2;
    tau3 = _tau3;
}

void OpenWristModel::set_positions(double _q1, double _q2, double _q3) {
    q1 = _q1;
    q2 = _q2;
    q3 = _q3;
    q1d_q1 = Integrator(q1);
    q2d_q2 = Integrator(q2);
    q3d_q3 = Integrator(q3);
}

void OpenWristModel::set_velocities(double _q1d, double _q2d, double _q3d) {
    q1d = _q1d;
    q2d = _q2d;
    q3d = _q3d;
    q1dd_q1d = Integrator(q1d);
    q2dd_q2d = Integrator(q2d);
    q3dd_q3d = Integrator(q3d);
}

void OpenWristModel::reset() {
    set_torques(0,0,0);
    set_positions(0,0,0);
    set_velocities(0,0,0);
}