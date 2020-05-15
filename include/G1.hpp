#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_G1(const std::vector<double>& qs){
	double qf = qs[0];
	double l1 = qs[1];
	double l2 = qs[2];
	double l3 = qs[3];
	double theta1 = qs[4];
	double theta2 = qs[5];
	double theta3 = qs[6];
	double sinqf = sin(qf);
	double cosqf = cos(qf);
	double sintheta1 = sin(theta1);
	double costheta1 = cos(theta1);
	double sintheta2 = sin(theta2);
	double costheta2 = cos(theta2);
	double sintheta3 = sin(theta3);
	double costheta3 = cos(theta3);

	double G1 = cos(qf+9.451666499990097E-2)*costheta2*1.108465732158948E-1-cos(qf+9.451666499990097E-2)*costheta3*1.108465732158948E-1-sin(qf+9.451666499990097E-2)*costheta1*1.279945977694297E-1+sin(qf+9.451666499990097E-2)*costheta2*6.399729888471484E-2+sin(qf+9.451666499990097E-2)*costheta3*6.399729888471484E-2+l3*costheta3*(cosqf*8.149720060173422E-1-sinqf*5.795003273524344E-1)*1.193096483708359-l2*cos(qf+9.451666499990097E-2)*costheta2*1.259064444922842+l3*cos(qf+9.451666499990097E-2)*costheta3*1.259064444922842+l1*sin(qf+9.451666499990097E-2)*costheta1*1.453842392395018-l2*sin(qf+9.451666499990097E-2)*costheta2*7.269211961975088E-1-l3*sin(qf+9.451666499990097E-2)*costheta3*7.269211961975088E-1-l2*costheta2*(cosqf*9.093480079900473E-1+sinqf*4.160362969123526E-1)*1.193096483708359+l1*costheta1*(cosqf*9.437600198180007E-2+sinqf*9.95536624264787E-1)*1.193096483708359;
	return G1;
}