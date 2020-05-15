#include "MeiiModel.hpp"
#include "Eqns.hpp"

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector4d;
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
    lim4(tau3_mot_cont * eta3, tau3_mot_max * eta3, seconds(2)),
    p(12),
    V(5,5),
    M(5,5),
    G(5),
    Tau(5),
    A(5,5),
    b(5),
    x(5)
{
    reset();
}

void MeiiModel::update(Time t)
{
    // get list of qs and q_dots
    std::vector<double> qs = {q1, q2, q3, q4, q5, q6, q7, q8, q1d, q2d, q3d, q4d, q5d, q6d, q7d, q8d};

    VectorXd qdot(5);
    qdot(0) = q1d;
    qdot(1) = q2d;
    qdot(2) = q3d;
    qdot(3) = q4d;
    qdot(4) = q5d;

    // Coriolis vector
    Clock MatCalcClock;
    if(threadpool){

        // V Matrix
        Clock SetupTime;
        v00 = p.push([&qs](int){return get_V11(qs);});
        v01 = p.push([&qs](int){return get_V12(qs);});
        v02 = p.push([&qs](int){return get_V13(qs);});
        v03 = p.push([&qs](int){return get_V14(qs);});
        v04 = p.push([&qs](int){return get_V15(qs);});
        v10 = p.push([&qs](int){return get_V21(qs);});
        v11 = p.push([&qs](int){return get_V22(qs);});
        v12 = p.push([&qs](int){return get_V23(qs);});
        v13 = p.push([&qs](int){return get_V24(qs);});
        v14 = p.push([&qs](int){return get_V25(qs);});
        v20 = p.push([&qs](int){return get_V31(qs);});
        v21 = p.push([&qs](int){return get_V32(qs);});
        v22 = p.push([&qs](int){return get_V33(qs);});
        v23 = p.push([&qs](int){return get_V34(qs);});
        v24 = p.push([&qs](int){return get_V35(qs);});
        v30 = p.push([&qs](int){return get_V41(qs);});
        v31 = p.push([&qs](int){return get_V42(qs);});
        v32 = p.push([&qs](int){return get_V43(qs);});
        v33 = p.push([&qs](int){return get_V44(qs);});
        v34 = p.push([&qs](int){return get_V45(qs);});
        v40 = p.push([&qs](int){return get_V51(qs);});
        v41 = p.push([&qs](int){return get_V52(qs);});
        v42 = p.push([&qs](int){return get_V53(qs);});
        v43 = p.push([&qs](int){return get_V54(qs);});
        v44 = p.push([&qs](int){return get_V55(qs);});

        // M Matrix
        m00 = p.push([&qs](int){return get_M11(qs);});
        m01 = p.push([&qs](int){return get_M12(qs);});
        m02 = p.push([&qs](int){return get_M13(qs);});
        m03 = p.push([&qs](int){return get_M14(qs);});
        m04 = p.push([&qs](int){return get_M15(qs);});
        m10 = p.push([&qs](int){return get_M21(qs);});
        m11 = p.push([&qs](int){return get_M22(qs);});
        m12 = p.push([&qs](int){return get_M23(qs);});
        m13 = p.push([&qs](int){return get_M24(qs);});
        m14 = p.push([&qs](int){return get_M25(qs);});
        m20 = p.push([&qs](int){return get_M31(qs);});
        m21 = p.push([&qs](int){return get_M32(qs);});
        m22 = p.push([&qs](int){return get_M33(qs);});
        m23 = p.push([&qs](int){return get_M34(qs);});
        m24 = p.push([&qs](int){return get_M35(qs);});
        m30 = p.push([&qs](int){return get_M41(qs);});
        m31 = p.push([&qs](int){return get_M42(qs);});
        m32 = p.push([&qs](int){return get_M43(qs);});
        m33 = p.push([&qs](int){return get_M44(qs);});
        m34 = p.push([&qs](int){return get_M45(qs);});
        m40 = p.push([&qs](int){return get_M51(qs);});
        m41 = p.push([&qs](int){return get_M52(qs);});
        m42 = p.push([&qs](int){return get_M53(qs);});
        m43 = p.push([&qs](int){return get_M54(qs);});
        m44 = p.push([&qs](int){return get_M55(qs);});

        // G Vector
        g0 = p.push([&qs](int){return get_G1(qs);});
        g1 = p.push([&qs](int){return get_G2(qs);});
        g2 = p.push([&qs](int){return get_G3(qs);});
        g3 = p.push([&qs](int){return get_G4(qs);});
        g4 = p.push([&qs](int){return get_G5(qs);});

        setup_time = double(SetupTime.get_elapsed_time().as_microseconds());
        Clock CompTime;

        Tau[0] = tau1 + hardstop_torque(q1,q1d,q1min,q1max,Khard1,Bhard1);
        Tau[1] = tau2 + hardstop_torque(q2,q2d,q2min,q2max,Khard1,Bhard1);
        Tau[2] = tau3 + hardstop_torque(q3,q3d,q3min,q3max,Khard,Bhard);
        Tau[3] = tau4 + hardstop_torque(q4,q4d,q4min,q4max,Khard,Bhard);
        Tau[4] = tau5 + hardstop_torque(q5,q5d,q5min,q5max,Khard,Bhard);

        G(0) = g0.get();
        G(1) = g1.get();
        G(2) = g2.get();
        G(3) = g3.get();
        G(4) = g4.get();
        
        M(0,0) = m00.get();
        M(0,1) = m01.get();
        M(0,2) = m02.get();
        M(0,3) = m03.get();
        M(0,4) = m04.get();
        M(1,0) = m10.get();
        M(1,1) = m11.get();
        M(1,2) = m12.get();
        M(1,3) = m13.get();
        M(1,4) = m14.get();
        M(2,0) = m20.get();
        M(2,1) = m21.get();
        M(2,2) = m22.get();
        M(2,3) = m23.get();
        M(2,4) = m24.get();
        M(3,0) = m30.get();
        M(3,1) = m31.get();
        M(3,2) = m32.get();
        M(3,3) = m33.get();
        M(3,4) = m34.get();
        M(4,0) = m40.get();
        M(4,1) = m41.get();
        M(4,2) = m42.get();
        M(4,3) = m43.get();
        M(4,4) = m44.get();

        // std::cout << "get_V" << std::endl;
        V(0,0) = v00.get();
        V(0,1) = v01.get();
        V(0,2) = v02.get();
        V(0,3) = v03.get();
        V(0,4) = v04.get();
        V(1,0) = v10.get();
        V(1,1) = v11.get();
        V(1,2) = v12.get();
        V(1,3) = v13.get();
        V(1,4) = v14.get();
        V(2,0) = v20.get();
        V(2,1) = v21.get();
        V(2,2) = v22.get();
        V(2,3) = v23.get();
        V(2,4) = v24.get();
        V(3,0) = v30.get();
        V(3,1) = v31.get();
        V(3,2) = v32.get();
        V(3,3) = v33.get();
        V(3,4) = v34.get();
        V(4,0) = v40.get();
        V(4,1) = v41.get();
        V(4,2) = v42.get();
        V(4,3) = v43.get();
        V(4,4) = v44.get();

        comp_time = double(CompTime.get_elapsed_time().as_microseconds());
    }
    else{
        V(0,0) = get_V11(qs);
        V(0,1) = get_V12(qs);
        V(0,2) = get_V13(qs);
        V(0,3) = get_V14(qs);
        V(0,4) = get_V15(qs);
        V(1,0) = get_V21(qs);
        V(1,1) = get_V22(qs);
        V(1,2) = get_V23(qs);
        V(1,3) = get_V24(qs);
        V(1,4) = get_V25(qs);
        V(2,0) = get_V31(qs);
        V(2,1) = get_V32(qs);
        V(2,2) = get_V33(qs);
        V(2,3) = get_V34(qs);
        V(2,4) = get_V35(qs);
        V(3,0) = get_V41(qs);
        V(3,1) = get_V42(qs);
        V(3,2) = get_V43(qs);
        V(3,3) = get_V44(qs);
        V(3,4) = get_V45(qs);
        V(4,0) = get_V51(qs);
        V(4,1) = get_V52(qs);
        V(4,2) = get_V53(qs);
        V(4,3) = get_V54(qs);
        V(4,4) = get_V55(qs);

        M(0,0) = get_M11(qs);
        M(0,1) = get_M12(qs);
        M(0,2) = get_M13(qs);
        M(0,3) = get_M14(qs);
        M(0,4) = get_M15(qs);
        M(1,0) = get_M21(qs);
        M(1,1) = get_M22(qs);
        M(1,2) = get_M23(qs);
        M(1,3) = get_M24(qs);
        M(1,4) = get_M25(qs);
        M(2,0) = get_M31(qs);
        M(2,1) = get_M32(qs);
        M(2,2) = get_M33(qs);
        M(2,3) = get_M34(qs);
        M(2,4) = get_M35(qs);
        M(3,0) = get_M41(qs);
        M(3,1) = get_M42(qs);
        M(3,2) = get_M43(qs);
        M(3,3) = get_M44(qs);
        M(3,4) = get_M45(qs);
        M(4,0) = get_M51(qs);
        M(4,1) = get_M52(qs);
        M(4,2) = get_M53(qs);
        M(4,3) = get_M54(qs);
        M(4,4) = get_M55(qs);

        G[0] = get_G1(qs);
        G[1] = get_G2(qs);
        G[2] = get_G3(qs);
        G[3] = get_G4(qs);
        G[4] = get_G5(qs);

        Tau[0] = tau1 + hardstop_torque(q1,q1d,q1min,q1max,Khard1,Bhard1);
        Tau[1] = tau2 + hardstop_torque(q2,q2d,q2min,q2max,Khard1,Bhard1);
        Tau[2] = tau3 + hardstop_torque(q3,q3d,q3min,q3max,Khard,Bhard);
        Tau[3] = tau4 + hardstop_torque(q4,q4d,q4min,q4max,Khard,Bhard);
        Tau[4] = tau5 + hardstop_torque(q5,q5d,q5min,q5max,Khard,Bhard);

        setup_time = 0.0;
        comp_time = 0.0;
    }

    // std::cout << "M:\n";
    // std::cout << M(0,0) << ", " << M(0,1) << ", " << M(0,2) << ", " << M(0,3) << ", " << M(0,4) << std::endl;
    // std::cout << M(1,0) << ", " << M(1,1) << ", " << M(1,2) << ", " << M(1,3) << ", " << M(1,4) << std::endl;
    // std::cout << M(2,0) << ", " << M(2,1) << ", " << M(2,2) << ", " << M(2,3) << ", " << M(2,4) << std::endl;
    // std::cout << M(3,0) << ", " << M(3,1) << ", " << M(3,2) << ", " << M(3,3) << ", " << M(3,4) << std::endl;
    // std::cout << M(4,0) << ", " << M(4,1) << ", " << M(4,2) << ", " << M(4,3) << ", " << M(4,4) << std::endl;

    // std::cout << "V:\n";
    // std::cout << V(0,0) << ", " << V(0,1) << ", " << V(0,2) << ", " << V(0,3) << ", " << V(0,4) << std::endl;
    // std::cout << V(1,0) << ", " << V(1,1) << ", " << V(1,2) << ", " << V(1,3) << ", " << V(1,4) << std::endl;
    // std::cout << V(2,0) << ", " << V(2,1) << ", " << V(2,2) << ", " << V(2,3) << ", " << V(2,4) << std::endl;
    // std::cout << V(3,0) << ", " << V(3,1) << ", " << V(3,2) << ", " << V(3,3) << ", " << V(3,4) << std::endl;
    // std::cout << V(4,0) << ", " << V(4,1) << ", " << V(4,2) << ", " << V(4,3) << ", " << V(4,4) << std::endl;

    // std::cout << "G:\n";
    // std::cout << G(0) << std::endl;
    // std::cout << G(1) << std::endl;
    // std::cout << G(2) << std::endl;
    // std::cout << G(3) << std::endl;
    // std::cout << G(4) << std::endl;

    mat_calc_time = MatCalcClock.get_elapsed_time().as_microseconds();

    A = M;// + M_mot;
    b = Tau - V*qdot - G;// - B - Fk;
    x = A.householderQr().solve(b);

    q1dd = x[0];
    q2dd = x[1];
    q3dd = x[2];
    q4dd = x[3];
    q5dd = x[4];

    // integrate acclerations to find velocities
    q1d = q1dd_q1d.update(q1dd, t);
    q2d = q2dd_q2d.update(q2dd, t);
    q3d = q3dd_q3d.update(q3dd, t);
    q4d = q4dd_q4d.update(q4dd, t);
    q5d = q5dd_q5d.update(q5dd, t);

    // integrate velocities to find positions
    q1 = q1d_q1.update(q1d, t);
    q2 = q2d_q2.update(q2d, t);
    q3 = q3d_q3.update(q3d, t);
    q4 = q4d_q4.update(q4d, t);
    q5 = q5d_q5.update(q5d, t);

    calc_dependent_joint_values();
}

void MeiiModel::set_torques(double _tau1, double _tau2, double _tau3, double _tau4, double _tau5) {
    tau1 = _tau1;
    tau2 = _tau2;
    tau3 = _tau3;
    tau4 = _tau4;
    tau5 = _tau5;
}

void MeiiModel::set_positions(double _q1, double _q2, double _q3, double _q4, double _q5, double _q6, double _q7, double _q8) {
    q1 = _q1;
    q2 = _q2;
    q3 = _q3;
    q4 = _q4;
    q5 = _q5;
    q1d_q1 = Integrator(q1);
    q2d_q2 = Integrator(q2);
    q3d_q3 = Integrator(q3);
    q4d_q4 = Integrator(q4);
    q5d_q5 = Integrator(q5);

    q6 = _q6;
    q7 = _q7;
    q8 = _q8;
}

void MeiiModel::set_velocities(double _q1d, double _q2d, double _q3d, double _q4d, double _q5d, double _q6d, double _q7d, double _q8d) {
    q1d = _q1d;
    q2d = _q2d;
    q3d = _q3d;
    q4d = _q4d;
    q5d = _q5d;
    q1dd_q1d = Integrator(q1d);
    q2dd_q2d = Integrator(q2d);
    q3dd_q3d = Integrator(q3d);
    q4dd_q4d = Integrator(q4d);
    q5dd_q5d = Integrator(q5d);

    q6d = _q6d;
    q7d = _q7d;
    q8d = _q8d;
}

void MeiiModel::reset() {
    set_torques(0,0,0,0,0);
    set_positions(0.0,0.0,0.1,0.1,0.1,1.0284436869069115694230731605785,1.0284436869069115694230731605785,1.0284436869069115694230731605785);
    set_velocities(0,0,0,0,0,0,0,0);
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
        jac(0,0) = (q3*(3*R*sin(t1)-(3*q4*sin(t1+t2))/2+(q4*sin(t1-t2))/2+sqrt(3)*a56*sin(t1)))/(2*sqrt(3*R*R+3*a56*a56+q3*q3+q4*q4-(q3*q4*cos(t1-t2))/2+(3*q3*q4*cos(t1+t2))/2-3*R*q3*cos(t1)-3*R*q4*cos(t2)-sqrt(3)*a56*q3*cos(t1)+sqrt(3)*a56*q4*cos(t2)));
        jac(0,1) = -(q4*((3*q3*sin(t1+t2))/2-3*R*sin(t2)+(q3*sin(t1-t2))/2+sqrt(3)*a56*sin(t2)))/(2*sqrt(3*R*R+3*a56*a56+q3*q3+q4*q4-(q3*q4*cos(t1-t2))/2+(3*q3*q4*cos(t1+t2))/2-3*R*q3*cos(t1)-3*R*q4*cos(t2)-sqrt(3)*a56*q3*cos(t1)+sqrt(3)*a56*q4*cos(t2)));
        jac(0,2) = 0;
        jac(1,0) = 0;
        jac(1,1) = (q4*(3*R*sin(t2)-(3*q5*sin(t2+t3))/2+(q5*sin(t2-t3))/2+sqrt(3)*a56*sin(t2)))/(2*sqrt(3*R*R+3*a56*a56+q4*q4+q5*q5-(q4*q5*cos(t2-t3))/2+(3*q4*q5*cos(t2+t3))/2-3*R*q4*cos(t2)-3*R*q5*cos(t3)-sqrt(3)*a56*q4*cos(t2)+sqrt(3)*a56*q5*cos(t3)));
        jac(1,2) = -(q5*((3*q4*sin(t2+t3))/2-3*R*sin(t3)+(q4*sin(t2-t3))/2+sqrt(3)*a56*sin(t3)))/(2*sqrt(3*R*R+3*a56*a56+q4*q4+q5*q5-(q4*q5*cos(t2-t3))/2+(3*q4*q5*cos(t2+t3))/2-3*R*q4*cos(t2)-3*R*q5*cos(t3)-sqrt(3)*a56*q4*cos(t2)+sqrt(3)*a56*q5*cos(t3)));
        jac(2,0) = -(q3*((3*q5*sin(t1+t3))/2-3*R*sin(t1)-(q5*sin(t1-t3))/2+sqrt(3)*a56*sin(t1)))/(2*sqrt(3*R*R+3*a56*a56+q3*q3+q5*q5-(q3*q5*cos(t1-t3))/2+(3*q3*q5*cos(t1+t3))/2-3*R*q3*cos(t1)-3*R*q5*cos(t3)+sqrt(3)*a56*q3*cos(t1)-sqrt(3)*a56*q5*cos(t3)));
        jac(2,1) = 0;
        jac(2,2) = -(q5*((3*q3*sin(t1+t3))/2-3*R*sin(t3)+(q3*sin(t1-t3))/2-sqrt(3)*a56*sin(t3)))/(2*sqrt(3*R*R+3*a56*a56+q3*q3+q5*q5-(q3*q5*cos(t1-t3))/2+(3*q3*q5*cos(t1+t3))/2-3*R*q3*cos(t1)-3*R*q5*cos(t3)+sqrt(3)*a56*q3*cos(t1)-sqrt(3)*a56*q5*cos(t3)));
        eq[0] =  sqrt(3*R*R + 3*a56*a56 + q3*q3 + q4*q4 - (q3*q4*cos(t1 - t2))/2 + (3*q3*q4*cos(t1 + t2))/2 - 3*R*q3*cos(t1) - 3*R*q4*cos(t2) - sqrt(3)*a56*q3*cos(t1) + sqrt(3)*a56*q4*cos(t2)) - sqrt(3)*r;
        eq[1] = sqrt(3*R*R + 3*a56*a56 + q4*q4 + q5*q5 - (q4*q5*cos(t2 - t3))/2 + (3*q4*q5*cos(t2 + t3))/2 - 3*R*q4*cos(t2) - 3*R*q5*cos(t3) - sqrt(3)*a56*q4*cos(t2) + sqrt(3)*a56*q5*cos(t3)) - sqrt(3)*r;
        eq[2] = sqrt(3*R*R + 3*a56*a56 + q3*q3 + q5*q5 - (q3*q5*cos(t1 - t3))/2 + (3*q3*q5*cos(t1 + t3))/2 - 3*R*q3*cos(t1) - 3*R*q5*cos(t3) + sqrt(3)*a56*q3*cos(t1) - sqrt(3)*a56*q5*cos(t3)) - sqrt(3)*r;
        thetas = thetas - jac.householderQr().solve(eq);
        iter = iter + 1;
    }

    q6 = thetas[0]; q7 = thetas[1]; q8 = thetas[2];

    // Velocity calculation
    std::vector<double> qs = {q1, q2, q3, q4, q5, q6, q7, q8};
    VectorXd qds(5);
    qds(0) = q1d;
    qds(1) = q2d;
    qds(2) = q3d;
    qds(3) = q4d;
    qds(4) = q5d;

    MatrixXd rho = calculate_rho(qs);
    
    // std::cout << rho(0,0) << ", " << q6d << ", " << q7d << std::endl;

    VectorXd qd = rho*qds;

    q6d = qd[5];
    q7d = qd[6];
    q8d = qd[7];
    // std::cout << q5d << ", " << q6d << ", " << q7d << std::endl;

    // NEED TO CHANGE ALL OF THESE AHHHHHHHHHHHHH
    double px = q3*sin(q6)*(1.0/3.0)+q4*sin(q7)*(1.0/3.0)+q5*sin(q8)*(1.0/3.0);
    double py = R*cos(alpha_5)*(1.0/3.0)+a56*sin(alpha_5)*(1.0/3.0)-R*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-R*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-a56*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-a56*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-q3*cos(alpha_5)*cos(q6)*(1.0/3.0)+q4*cos(q7)*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)+q5*cos(q8)*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0);
    double pz = a56*cos(alpha_5)*(-1.0/3.0)+R*sin(alpha_5)*(1.0/3.0)-R*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)-R*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)+a56*(cos(alpha_5)*(1.0/2.0)-sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)+a56*(cos(alpha_5)*(1.0/2.0)+sqrt(3.0)*sin(alpha_5)*(1.0/2.0))*(1.0/3.0)-q3*sin(alpha_5)*cos(q6)*(1.0/3.0)+q4*cos(q7)*(sin(alpha_5)*(1.0/2.0)+sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0)+q5*cos(q8)*(sin(alpha_5)*(1.0/2.0)-sqrt(3.0)*cos(alpha_5)*(1.0/2.0))*(1.0/3.0);

    UnitVecs[0] = -(px*cos(alpha_13)*-3.0+sqrt(3.0)*px*sin(alpha_13)+q3*cos(alpha_13)*sin(q6)+q4*cos(alpha_13)*sin(q7)*2.0-sqrt(3.0)*q3*sin(alpha_13)*sin(q6))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[1] = (py*cos(alpha_13)*3.0-sqrt(3.0)*py*sin(alpha_13)+q3*cos(alpha_5)*cos(alpha_13)*cos(q6)-q4*cos(alpha_5)*cos(alpha_13)*cos(q7)+sqrt(3.0)*a56*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*R*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*R*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*a56*sin(alpha_5)*sin(alpha_13)-sqrt(3.0)*q3*cos(alpha_5)*sin(alpha_13)*cos(q6)+sqrt(3.0)*q4*cos(alpha_13)*sin(alpha_5)*cos(q7))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[2] = (pz*cos(alpha_13)*3.0-sqrt(3.0)*pz*sin(alpha_13)+q3*cos(alpha_13)*sin(alpha_5)*cos(q6)-q4*cos(alpha_13)*sin(alpha_5)*cos(q7)+sqrt(3.0)*R*cos(alpha_5)*cos(alpha_13)-sqrt(3.0)*a56*cos(alpha_5)*sin(alpha_13)+sqrt(3.0)*a56*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*R*sin(alpha_5)*sin(alpha_13)-sqrt(3.0)*q4*cos(alpha_5)*cos(alpha_13)*cos(q7)-sqrt(3.0)*q3*sin(alpha_5)*sin(alpha_13)*cos(q6))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[3] = (px*sin(alpha_13)*-3.0+q3*sin(alpha_13)*sin(q6)+q4*sin(alpha_13)*sin(q7)*2.0-sqrt(3.0)*px*cos(alpha_13)+sqrt(3.0)*q3*cos(alpha_13)*sin(q6))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[4] = -(py*sin(alpha_13)*3.0+sqrt(3.0)*py*cos(alpha_13)+q3*cos(alpha_5)*sin(alpha_13)*cos(q6)-q4*cos(alpha_5)*sin(alpha_13)*cos(q7)-sqrt(3.0)*R*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*a56*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*a56*cos(alpha_13)*sin(alpha_5)-sqrt(3.0)*R*sin(alpha_5)*sin(alpha_13)+sqrt(3.0)*q3*cos(alpha_5)*cos(alpha_13)*cos(q6)+sqrt(3.0)*q4*sin(alpha_5)*sin(alpha_13)*cos(q7))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));
    UnitVecs[5] = -(pz*sin(alpha_13)*3.0+sqrt(3.0)*pz*cos(alpha_13)+q3*sin(alpha_5)*sin(alpha_13)*cos(q6)-q4*sin(alpha_5)*sin(alpha_13)*cos(q7)+sqrt(3.0)*a56*cos(alpha_5)*cos(alpha_13)+sqrt(3.0)*R*cos(alpha_5)*sin(alpha_13)-sqrt(3.0)*R*cos(alpha_13)*sin(alpha_5)+sqrt(3.0)*a56*sin(alpha_5)*sin(alpha_13)+sqrt(3.0)*q3*cos(alpha_13)*sin(alpha_5)*cos(q6)-sqrt(3.0)*q4*cos(alpha_5)*sin(alpha_13)*cos(q7))/(r*(sqrt(3.0)*pow(cos(alpha_13),2.0)+sqrt(3.0)*pow(sin(alpha_13),2.0)));

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

    // std::cout << a1 << ", " << a2 << ", " << a3 << ", " << o1 << ", " << o2 << ", " << o3 << ", " << n1 << ", " << n2 << ", " << n3 << std::endl;

    double alpha = atan2(-n3,n1);
    double beta = atan2(n2,sqrt(n3*n3+n1*n1));
    double gamma = atan2(-a2,o2);

    q9  = px;
    q10  = py;
    q11 = pz;
    q12 = alpha;
    q13 = beta;
    q14 = gamma;

    // std::cout << q8 << ", " << q9 << ", " << q10 << ", " << q11 << ", " << q12 << ", " << q13 << std::endl;
}