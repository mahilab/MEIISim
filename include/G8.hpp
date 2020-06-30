#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_G8(const std::vector<double>& qs){
	double qe = qs[0];
	double qf = qs[1];
	double l1 = qs[2];
	double l2 = qs[3];
	double l3 = qs[4];
	double theta1 = qs[5];
	double theta2 = qs[6];
	double theta3 = qs[7];
	double sinqe = sin(qe);
	double cosqe = cos(qe);
	double sinqf = sin(qf);
	double cosqf = cos(qf);
	double sintheta1 = sin(theta1);
	double costheta1 = cos(theta1);
	double sintheta2 = sin(theta2);
	double costheta2 = cos(theta2);
	double sintheta3 = sin(theta3);
	double costheta3 = cos(theta3);

	double G8 = (l3-8.803883999999584E-2)*(costheta3*sinqe-cosqe*sintheta3*(cosqf*5.795003273562997E-1+sinqf*8.149720060185928E-1)*1.0)*1.45384239240002+l3*costheta3*sinqe*1.193096483699946-l3*cos(qf+9.451666500000044E-2)*cosqe*sintheta3*5.965482418499732E-1-l3*sin(qf+9.451666500000044E-2)*cosqe*sintheta3*1.033251864050044;
	return G8;
}