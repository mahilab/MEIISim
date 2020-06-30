#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_G3(const std::vector<double>& qs){
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

	double G3 = sinqe*sintheta1*2.646938876100194-cos(qf+9.451666500000044E-2)*cosqe*costheta1*1.193096483699946-cosqe*costheta1*(cosqf*9.955366242634227E-1-sinqf*9.437600198272378E-2)*1.45384239240002;
	return G3;
}