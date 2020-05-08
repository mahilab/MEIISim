#include "MeiiModel.hpp"
#include "Eqns.hpp"
#include <Eigen/Dense>
// #include <thread>

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::VectorXd;

using namespace mahi::util;

inline double hardstop_torque(double q, double qd, double qmin, double qmax, double K, double B) {
    if (q < qmin)
        return K * (qmin - q) - B * qd;
    else if (q > qmax){
        return K * (qmax - q) - B * qd;
    }
    else
        return 0;
}

MeiiModel::MeiiModel() : 
    lim1(tau1_mot_cont * eta1, tau1_mot_max * eta1, seconds(2)),
    lim2(tau2_mot_cont * eta2, tau2_mot_max * eta2, seconds(2)),
    lim3(tau3_mot_cont * eta3, tau3_mot_max * eta3, seconds(2)),
    p(12)
{
    reset();
    // std::cout << "here1";
    // std::cout << "here2";
}

void MeiiModel::update(Time t)
{
    // std::cout << "here3";
    // limit torques
    // tau1 = lim1.limit(tau1);
    // tau2 = lim2.limit(tau2);
    // tau3 = lim3.limit(tau3);

    // get list of qs and q_dots
    std::vector<double> qs = {q1, q2, q3, q4, q5, q6, q1d, q2d, q3d, q4d, q5d, q6d};

    Vector3d qdot(q1d,q2d,q3d);
    // Coriolis vector
    Matrix3d V;
    Matrix3d M;
    Vector3d G;
    Vector3d Tau;
    Clock MatCalcClock;
    if(threadpool){

        // V Matrix
        Clock SetupTime;

        std::vector<double> qs1 = qs;
        auto v00 = p.push([&qs1](int){return get_V11(qs1);});
        int padding1[32];

        std::vector<double> qs2 = qs;
        auto v01 = p.push([&qs2](int){return get_V12(qs2);});
        int padding2[32];

        std::vector<double> qs3 = qs;
        auto v02 = p.push([&qs3](int){return get_V13(qs3);});
        int padding3[32];

        std::vector<double> qs4 = qs;
        auto v10 = p.push([&qs4](int){return get_V21(qs4);});
        int padding4[32];

        std::vector<double> qs5 = qs;
        auto v11 = p.push([&qs5](int){return get_V22(qs5);});
        int padding5[32];

        std::vector<double> qs6 = qs;
        auto v12 = p.push([&qs6](int){return get_V23(qs6);});
        int padding6[32];

        std::vector<double> qs7 = qs;
        auto v20 = p.push([&qs7](int){return get_V31(qs7);});
        int padding7[32];

        std::vector<double> qs8 = qs;
        auto v21 = p.push([&qs8](int){return get_V32(qs8);});
        int padding8[32];

        std::vector<double> qs9 = qs;
        auto v22 = p.push([&qs9](int){return get_V33(qs9);});
        int padding9[32];

        // M Matrix

        std::vector<double> qs10 = qs;
        auto m00 = p.push([&qs10](int){return get_M11(qs10);});
        int padding10[32];

        std::vector<double> qs11 = qs;
        auto m01 = p.push([&qs11](int){return get_M12(qs11);});
        int padding11[32];

        std::vector<double> qs12 = qs;
        auto m02 = p.push([&qs12](int){return get_M13(qs12);});
        int padding12[32];

        std::vector<double> qs13 = qs;
        auto m10 = p.push([&qs13](int){return get_M21(qs13);});
        int padding13[32];

        std::vector<double> qs14 = qs;
        auto m11 = p.push([&qs14](int){return get_M22(qs14);});
        int padding14[32];

        std::vector<double> qs15 = qs;
        auto m12 = p.push([&qs15](int){return get_M23(qs15);});
        int padding15[32];

        std::vector<double> qs16 = qs;
        auto m20 = p.push([&qs16](int){return get_M31(qs16);});
        int padding16[32];

        std::vector<double> qs17 = qs;
        auto m21 = p.push([&qs17](int){return get_M32(qs17);});
        int padding17[32];

        std::vector<double> qs18 = qs;
        auto m22 = p.push([&qs18](int){return get_M33(qs18);});
        int padding18[32];

        // G Vector

        std::vector<double> qs19 = qs;
        auto g0 = p.push([&qs19](int){return get_G1(qs19);});
        int padding19[32];

        std::vector<double> qs20 = qs;
        auto g1 = p.push([&qs20](int){return get_G2(qs20);});
        int padding20[32];

        std::vector<double> qs21 = qs;
        auto g2 = p.push([&qs21](int){return get_G3(qs21);});
        int padding21[32];

        setup_time = double(SetupTime.get_elapsed_time().as_microseconds());
        Clock CompTime;

        Tau[0] = tau1 + hardstop_torque(q1,q1d,q1min,q1max,Khard,Bhard);
        Tau[1] = tau2 + hardstop_torque(q2,q2d,q2min,q2max,Khard,Bhard);
        Tau[2] = tau3 + hardstop_torque(q3,q3d,q3min,q3max,Khard,Bhard);

        G(0) = g0.get();
        G(1) = g1.get();
        G(2) = g2.get();
        
        M(0,0) = m00.get();
        M(0,1) = m01.get();
        M(0,2) = m02.get();
        M(1,0) = m10.get();
        M(1,1) = m11.get();
        M(1,2) = m12.get();
        M(2,0) = m20.get();
        M(2,1) = m21.get();
        M(2,2) = m22.get();

        V(0,0) = v00.get();
        V(0,1) = v01.get();
        V(0,2) = v02.get();
        V(1,0) = v10.get();
        V(1,1) = v11.get();
        V(1,2) = v12.get();
        V(2,0) = v20.get();
        V(2,1) = v21.get();
        V(2,2) = v22.get();

        comp_time = double(CompTime.get_elapsed_time().as_microseconds());
    }
    else{
        V(0,0) = get_V11(qs);
        V(0,1) = get_V12(qs);
        V(0,2) = get_V13(qs);
        V(1,0) = get_V21(qs);
        V(1,1) = get_V22(qs);
        V(1,2) = get_V23(qs);
        V(2,0) = get_V31(qs);
        V(2,1) = get_V32(qs);
        V(2,2) = get_V33(qs);
        
        M(0,0) = get_M11(qs);
        M(0,1) = get_M12(qs);
        M(0,2) = get_M13(qs);
        M(1,0) = get_M21(qs);
        M(1,1) = get_M22(qs);
        M(1,2) = get_M23(qs);
        M(2,0) = get_M31(qs);
        M(2,1) = get_M32(qs);
        M(2,2) = get_M33(qs);

        G[0] = get_G1(qs);
        G[1] = get_G2(qs);
        G[2] = get_G3(qs);

        Tau[0] = tau1 + hardstop_torque(q1,q1d,q1min,q1max,Khard,Bhard);
        Tau[1] = tau2 + hardstop_torque(q2,q2d,q2min,q2max,Khard,Bhard);
        Tau[2] = tau3 + hardstop_torque(q3,q3d,q3min,q3max,Khard,Bhard);

        setup_time = 0.0;
        comp_time = 0.0;
    }

    mat_calc_time = MatCalcClock.get_elapsed_time().as_microseconds();

    Matrix3d A = M;// + M_mot;
    Vector3d b = Tau - V*qdot - G;// - B - Fk;
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

    calc_dependent_joint_values();
}

void MeiiModel::set_torques(double _tau1, double _tau2, double _tau3) {
    tau1 = _tau1;
    tau2 = _tau2;
    tau3 = _tau3;
}

void MeiiModel::set_positions(double _q1, double _q2, double _q3, double _q4, double _q5, double _q6) {
    q1 = _q1;
    q2 = _q2;
    q3 = _q3;
    q1d_q1 = Integrator(q1);
    q2d_q2 = Integrator(q2);
    q3d_q3 = Integrator(q3);

    q4 = _q4;
    q5 = _q5;
    q6 = _q6;
}

void MeiiModel::set_velocities(double _q1d, double _q2d, double _q3d, double _q4d, double _q5d, double _q6d) {
    q1d = _q1d;
    q2d = _q2d;
    q3d = _q3d;
    q1dd_q1d = Integrator(q1d);
    q2dd_q2d = Integrator(q2d);
    q3dd_q3d = Integrator(q3d);

    q4d = _q4d;
    q5d = _q5d;
    q6d = _q6d;
}

void MeiiModel::reset() {
    set_torques(0,0,0);
    set_positions(0.1,0.1,0.1,1.0284436869069115694230731605785,1.0284436869069115694230731605785,1.0284436869069115694230731605785);
    set_velocities(0,0,0,0,0,0);
}

void MeiiModel::calc_dependent_joint_values() {
    Vector3d thetas(PI/4, PI/4, PI/4);
    Vector3d old_thetas(PI, PI, PI);
    size_t iter = 0;
    Matrix3d jac;
    Vector3d eq;
    VectorXd UnitVecs(6);

    while (rms(old_thetas-thetas) > 1e-12 && iter < 20){
        old_thetas = thetas;
        double t1 = thetas[0]; double t2 = thetas[1]; double t3 = thetas[2];
        jac(0,0) = (q1*(3*R*sin(t1)-(3*q2*sin(t1+t2))/2+(q2*sin(t1-t2))/2+sqrt(3)*a56*sin(t1)))/(2*sqrt(3*R*R+3*a56*a56+q1*q1+q2*q2-(q1*q2*cos(t1-t2))/2+(3*q1*q2*cos(t1+t2))/2-3*R*q1*cos(t1)-3*R*q2*cos(t2)-sqrt(3)*a56*q1*cos(t1)+sqrt(3)*a56*q2*cos(t2)));
        jac(0,1) = -(q2*((3*q1*sin(t1+t2))/2-3*R*sin(t2)+(q1*sin(t1-t2))/2+sqrt(3)*a56*sin(t2)))/(2*sqrt(3*R*R+3*a56*a56+q1*q1+q2*q2-(q1*q2*cos(t1-t2))/2+(3*q1*q2*cos(t1+t2))/2-3*R*q1*cos(t1)-3*R*q2*cos(t2)-sqrt(3)*a56*q1*cos(t1)+sqrt(3)*a56*q2*cos(t2)));
        jac(0,2) = 0;
        jac(1,0) = 0;
        jac(1,1) = (q2*(3*R*sin(t2)-(3*q3*sin(t2+t3))/2+(q3*sin(t2-t3))/2+sqrt(3)*a56*sin(t2)))/(2*sqrt(3*R*R+3*a56*a56+q2*q2+q3*q3-(q2*q3*cos(t2-t3))/2+(3*q2*q3*cos(t2+t3))/2-3*R*q2*cos(t2)-3*R*q3*cos(t3)-sqrt(3)*a56*q2*cos(t2)+sqrt(3)*a56*q3*cos(t3)));
        jac(1,2) = -(q3*((3*q2*sin(t2+t3))/2-3*R*sin(t3)+(q2*sin(t2-t3))/2+sqrt(3)*a56*sin(t3)))/(2*sqrt(3*R*R+3*a56*a56+q2*q2+q3*q3-(q2*q3*cos(t2-t3))/2+(3*q2*q3*cos(t2+t3))/2-3*R*q2*cos(t2)-3*R*q3*cos(t3)-sqrt(3)*a56*q2*cos(t2)+sqrt(3)*a56*q3*cos(t3)));
        jac(2,0) = -(q1*((3*q3*sin(t1+t3))/2-3*R*sin(t1)-(q3*sin(t1-t3))/2+sqrt(3)*a56*sin(t1)))/(2*sqrt(3*R*R+3*a56*a56+q1*q1+q3*q3-(q1*q3*cos(t1-t3))/2+(3*q1*q3*cos(t1+t3))/2-3*R*q1*cos(t1)-3*R*q3*cos(t3)+sqrt(3)*a56*q1*cos(t1)-sqrt(3)*a56*q3*cos(t3)));
        jac(2,1) = 0;
        jac(2,2) = -(q3*((3*q1*sin(t1+t3))/2-3*R*sin(t3)+(q1*sin(t1-t3))/2-sqrt(3)*a56*sin(t3)))/(2*sqrt(3*R*R+3*a56*a56+q1*q1+q3*q3-(q1*q3*cos(t1-t3))/2+(3*q1*q3*cos(t1+t3))/2-3*R*q1*cos(t1)-3*R*q3*cos(t3)+sqrt(3)*a56*q1*cos(t1)-sqrt(3)*a56*q3*cos(t3)));
        eq[0] =  sqrt(3*R*R + 3*a56*a56 + q1*q1 + q2*q2 - (q1*q2*cos(t1 - t2))/2 + (3*q1*q2*cos(t1 + t2))/2 - 3*R*q1*cos(t1) - 3*R*q2*cos(t2) - sqrt(3)*a56*q1*cos(t1) + sqrt(3)*a56*q2*cos(t2)) - sqrt(3)*r;
        eq[1] = sqrt(3*R*R + 3*a56*a56 + q2*q2 + q3*q3 - (q2*q3*cos(t2 - t3))/2 + (3*q2*q3*cos(t2 + t3))/2 - 3*R*q2*cos(t2) - 3*R*q3*cos(t3) - sqrt(3)*a56*q2*cos(t2) + sqrt(3)*a56*q3*cos(t3)) - sqrt(3)*r;
        eq[2] = sqrt(3*R*R + 3*a56*a56 + q1*q1 + q3*q3 - (q1*q3*cos(t1 - t3))/2 + (3*q1*q3*cos(t1 + t3))/2 - 3*R*q1*cos(t1) - 3*R*q3*cos(t3) + sqrt(3)*a56*q1*cos(t1) - sqrt(3)*a56*q3*cos(t3)) - sqrt(3)*r;
        thetas = thetas - jac.householderQr().solve(eq);
        iter = iter + 1;
    }

    q4 = thetas[0]; q5 = thetas[1]; q6 = thetas[2];

    // Velocity calculation
    std::vector<double> qs = {q1, q2, q3, q4, q5, q6};
    Vector3d qds(q1d, q2d, q3d);

    MatrixXd rho = calculate_rho(qs);
    
    VectorXd qd = rho*qds;

    q4d = qd[3];
    q5d = qd[4];
    q6d = qd[5];

    double px = q1*sin(q4)*(1.0/3.0)+q2*sin(q5)*(1.0/3.0)+q3*sin(q6)*(1.0/3.0);
    double py = R*cos(alpha_5)*(1.0/3.0)+a56*sin(alpha_5)*(1.0/3.0)-R*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-R*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-a56*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-a56*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-q1*cos(alpha_5)*cos(q4)*(1.0/3.0)+q2*cos(q5)*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)+q3*cos(q6)*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0);
    double pz = a56*cos(alpha_5)*(-1.0/3.0)+R*sin(alpha_5)*(1.0/3.0)-R*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-R*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)+a56*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)+a56*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-q1*sin(alpha_5)*cos(q4)*(1.0/3.0)+q2*cos(q5)*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)+q3*cos(q6)*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0);

    UnitVecs[0] = -(px*cos(alpha_13)*-3.0+sqrt(3.0)*px*sin(alpha_13)+q1*cos(alpha_13)*sin(q4)+q2*cos(alpha_13)*sin(q5)*2.0-sqrt(3.0)*q1*sin(alpha_13)*sin(q4))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[1] = (py*cos(alpha_13)*3.0-sqrt(3.0)*py*sin(alpha_13)+q1*cos(alpha_5)*cos(alpha_13)*cos(q4)-q2*cos(alpha_5)*cos(alpha_13)*cos(q5)+sqrt(3.0)*a56*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*R*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*R*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*a56*sin(alpha_5)*sin(alpha_13)-sqrt(3.0)*q1*cos(alpha_5)*sin(alpha_13)*cos(q4)+sqrt(3.0)*q2*cos(alpha_13)*sin(alpha_5)*cos(q5))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[2] = (pz*cos(alpha_13)*3.0-sqrt(3.0)*pz*sin(alpha_13)+q1*cos(alpha_13)*sin(alpha_5)*cos(q4)-q2*cos(alpha_13)*sin(alpha_5)*cos(q5)+sqrt(3.0)*R*cos(alpha_5)*cos(alpha_13)-sqrt(3.0)*a56*cos(alpha_5)*sin(alpha_13)+sqrt(3.0)*a56*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*R*sin(alpha_5)*sin(alpha_13)-sqrt(3.0)*q2*cos(alpha_5)*cos(alpha_13)*cos(q5)-sqrt(3.0)*q1*sin(alpha_5)*sin(alpha_13)*cos(q4))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[3] = (px*sin(alpha_13)*-3.0+q1*sin(alpha_13)*sin(q4)+q2*sin(alpha_13)*sin(q5)*2.0-sqrt(3.0)*px*cos(alpha_13)+sqrt(3.0)*q1*cos(alpha_13)*sin(q4))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[4] = -(py*sin(alpha_13)*3.0+sqrt(3.0)*py*cos(alpha_13)+q1*cos(alpha_5)*sin(alpha_13)*cos(q4)-q2*cos(alpha_5)*sin(alpha_13)*cos(q5)-sqrt(3.0)*R*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*a56*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*a56*cos(alpha_13)*sin(alpha_5)-sqrt(3.0)*R*sin(alpha_5)*sin(alpha_13)+sqrt(3.0)*q1*cos(alpha_5)*cos(alpha_13)*cos(q4)+sqrt(3.0)*q2*sin(alpha_5)*sin(alpha_13)*cos(q5))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[5] = -(pz*sin(alpha_13)*3.0+sqrt(3.0)*pz*cos(alpha_13)+q1*sin(alpha_5)*sin(alpha_13)*cos(q4)-q2*sin(alpha_5)*sin(alpha_13)*cos(q5)+sqrt(3.0)*a56*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*R*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*R*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*a56*sin(alpha_5)*sin(alpha_13)+sqrt(3.0)*q1*cos(alpha_13)*sin(alpha_5)*cos(q4)-sqrt(3.0)*q2*cos(alpha_5)*sin(alpha_13)*cos(q5))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));

    UnitVecs.head(3) = UnitVecs.head(3)/UnitVecs.head(3).norm();
    UnitVecs.tail<3>() = UnitVecs.tail<3>()/UnitVecs.tail<3>().norm();

    double a1 = UnitVecs(0);
    double a2 = UnitVecs(1);
    double a3 = UnitVecs(2);
    double o1 = UnitVecs(3);
    double o2 = UnitVecs(4);
    double o3 = UnitVecs(5);
    double n1 =  o2*a3 - a2*o3;
    double n2 = -o1*a3 + a1*o3;
    double n3 =  o1*a2 - o2*a1;

    double alpha = atan2(-n3,n1);
    double beta = atan2(n2,sqrt(n3*n3+n1*n1));
    double gamma = atan2(-a2,o2);

    q7  = px;
    q8  = py;
    q9  = pz;
    q10 = alpha;
    q11 = beta;
    q12 = gamma;
}