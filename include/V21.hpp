#include <Mahi/Util.hpp>
using namespace mahi::util;

double get_V21(std::vector<double> qs){
	double l1 = qs[0];
	double l2 = qs[1];
	double l3 = qs[2];
	double theta1 = qs[3];
	double theta2 = qs[4];
	double theta3 = qs[5];
	double l1_dot = qs[6];
	double l2_dot = qs[7];
	double l3_dot = qs[8];
	double theta1_dot = qs[9];
	double theta2_dot = qs[10];
	double theta3_dot = qs[11];
	double sintheta1 = sin(theta1);
	double costheta1 = cos(theta1);
	double sintheta2 = sin(theta2);
	double costheta2 = cos(theta2);
	double sintheta3 = sin(theta3);
	double costheta3 = cos(theta3);

	return V21;
}