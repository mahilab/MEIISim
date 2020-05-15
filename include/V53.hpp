#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_V53(const std::vector<double>& qs){
	double qe = qs[0];
	double qf = qs[1];
	double l1 = qs[2];
	double l2 = qs[3];
	double l3 = qs[4];
	double theta1 = qs[5];
	double theta2 = qs[6];
	double theta3 = qs[7];
	double qe_dot = qs[8];
	double qf_dot = qs[9];
	double l1_dot = qs[10];
	double l2_dot = qs[11];
	double l3_dot = qs[12];
	double theta1_dot = qs[13];
	double theta2_dot = qs[14];
	double theta3_dot = qs[15];
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

	return V53;
}