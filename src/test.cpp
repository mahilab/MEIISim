#include <Mahi/Util.hpp>
#include <MeiiModel.hpp>
// #include <Eqns.hpp>
#include <Eigen/Dense>

using namespace mahi::util;
using Eigen::Matrix3d;
using Eigen::Vector3d;

int main(){
    // std::vector<double> qs = {0.1, 0.1, 0.11, 1.0390, 0.9300, 1.0663, 0.001, 0.0001, -0.0003, 0.001, 0.001, -0.003};
    std::vector<double> qs = {0.11, 0.1, 0.09, 1.0284, 1.0284, 1.0284, 0.002, 0.002, 0.002, 0.000, 0.000, 0.000};

    Clock functime;
    // Mass Matrix
    // Matrix3d M;
    // M(0,0) = get_M11(qs);
    // M(0,1) = get_M12(qs);
    // M(0,2) = get_M13(qs);
    // M(1,0) = get_M21(qs);
    // M(1,1) = get_M22(qs);
    // M(1,2) = get_M23(qs);
    // M(2,0) = get_M31(qs);
    // M(2,1) = get_M32(qs);
    // M(2,2) = get_M33(qs);
    // // Coriolis vector
    // Matrix3d V;
    // V(0,0) = get_V11(qs);
    // V(0,1) = get_V12(qs);
    // V(0,2) = get_V13(qs);
    // V(1,0) = get_V21(qs);
    // V(1,1) = get_V22(qs);
    // V(1,2) = get_V23(qs);
    // V(2,0) = get_V31(qs);
    // V(2,1) = get_V32(qs);
    // V(2,2) = get_V33(qs);
    // // Gravity vector
    // Vector3d G;
    // G[0] = get_G1(qs);
    // G[1] = get_G2(qs);
    // G[2] = get_G3(qs);

    // Vector3d Tau;
    // Tau[0] = 0.0;
    // Tau[1] = 0.0;
    // Tau[2] = 0.0;

    // Vector3d Qdot;
    // Qdot[0] = qs[6];
    // Qdot[1] = qs[7];
    // Qdot[2] = qs[8];

    // Matrix3d A = M;// + M_mot;
    // Vector3d b = Tau - V*Qdot - G;// - B - Fk;
    // Vector3d x = A.householderQr().solve(b);
    // std::cout << functime.get_elapsed_time().as_microseconds() << std::endl;
    // // std::cout << M << std::endl;
    // std::cout << V*Qdot << std::endl;
    // std::cout << G << std::endl;
    // std::cout << M(0,0) << " " << M(0,1) << " " << M(0,2) << std::endl;
    // std::cout << M(1,0) << " " << M(1,1) << " " << M(1,2) << std::endl;
    // std::cout << M(2,0) << " " << M(2,1) << " " << M(2,2) << std::endl;
    // std::cout << x << std::endl;
    // std::cout << rms(G) << std::endl;
    // MeiiModel meii;
    // rho_test(0.1,0.1,0.1);

    return 0;
}