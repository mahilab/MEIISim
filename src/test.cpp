#include <Mahi/Util.hpp>
#include <Eqns.hpp>
using namespace mahi::util;

int main(){
    std::vector<double> qs = {0.1, 0.1, 0.11, 1.0390, 0.9300, 1.0663, 0.001, 0.0001, -0.0003, 0.001, 0.001, -0.003};

    // double qs1;
    // double qs2;
    // double qs3;
    // double qs4;
    // double qs5;
    // double qs6;
    // double qs7;
    // double qs8;
    // double qs9;
    // double qs10;
    // double qs11;
    // double qs12;

    // std::cin >> qs1;
    // std::cin >> qs2;
    // std::cin >> qs3;
    // std::cin >> qs4;
    // std::cin >> qs5;
    // std::cin >> qs6;
    // std::cin >> qs7;
    // std::cin >> qs8;
    // std::cin >> qs9;
    // std::cin >> qs10;
    // std::cin >> qs11;
    // std::cin >> qs12;

    // std::vector<double> qs = {qs1, qs2, qs3, qs4, qs5, qs6, qs7, qs8, qs9, qs10, qs11, qs12};

    Clock functime;
    double G1 = get_G1(qs);
    double G2 = get_G2(qs);
    double G3 = get_G3(qs);
    double M11 = get_M11(qs);
    double M12 = get_M12(qs);
    double M13 = get_M13(qs);
    double M21 = get_M21(qs);
    double M22 = get_M22(qs);
    double M23 = get_M23(qs);
    double M31 = get_M31(qs);
    double M32 = get_M32(qs);
    double M33 = get_M33(qs);
    double V11 = get_V11(qs);
    double V12 = get_V12(qs);
    double V13 = get_V13(qs);
    double V21 = get_V21(qs);
    double V22 = get_V22(qs);
    double V23 = get_V23(qs);
    double V31 = get_V31(qs);
    double V32 = get_V32(qs);
    double V33 = get_V33(qs);
    std::cout << functime.get_elapsed_time().as_microseconds() << std::endl;
    std::cout << G1 << std::endl;
    std::cout << G2 << std::endl;
    std::cout << G3 << std::endl;
    std::cout << M11 << std::endl;
    std::cout << M12 << std::endl;
    std::cout << M13 << std::endl;
    std::cout << M21 << std::endl;
    std::cout << M22 << std::endl;
    std::cout << M23 << std::endl;
    std::cout << M31 << std::endl;
    std::cout << M32 << std::endl;
    std::cout << M33 << std::endl;
    std::cout << V11 << std::endl;
    std::cout << V12 << std::endl;
    std::cout << V13 << std::endl;
    std::cout << V21 << std::endl;
    std::cout << V22 << std::endl;
    std::cout << V23 << std::endl;
    std::cout << V31 << std::endl;
    std::cout << V32 << std::endl;
    std::cout << V33 << std::endl;
    return 0;
}