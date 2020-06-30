#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_G2(const std::vector<double>& qs){
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

	double G2 = cosqe*(cosqf*5.795003273562997E-1+sinqf*8.149720060185928E-1)*-5.574031732461515E-4-cosqe*(cosqf*8.149720060185928E-1-sinqf*5.795003273562997E-1)*1.519201330992814E-1-cosqe*(cosqf*4.160362969070661E-1-sinqf*9.093480080013023E-1)*5.574031732461515E-4+cosqe*(cosqf*9.093480080013023E-1+sinqf*4.160362969070661E-1)*1.519201330992814E-1-cosqe*(cosqf*9.437600198272378E-2+sinqf*9.955366242634227E-1)*1.519201330992814E-1+cosqe*(cosqf*9.955366242634227E-1-sinqf*9.437600198272378E-2)*5.574031732461515E-4-l2*cos(qf+9.451666500000044E-2)*cosqe*costheta2*1.033251864050044+l3*cos(qf+9.451666500000044E-2)*cosqe*costheta3*1.033251864050044+l1*sin(qf+9.451666500000044E-2)*cosqe*costheta1*1.193096483699946-l2*sin(qf+9.451666500000044E-2)*cosqe*costheta2*5.965482418499732E-1-l3*sin(qf+9.451666500000044E-2)*cosqe*costheta3*5.965482418499732E-1+cosqe*costheta3*(l3-8.803883999999584E-2)*(cosqf*8.149720060185928E-1-sinqf*5.795003273562997E-1)*1.45384239240002-cosqe*costheta2*(l2-8.803883999999584E-2)*(cosqf*9.093480080013023E-1+sinqf*4.160362969070661E-1)*1.45384239240002+cosqe*costheta1*(l1-8.803883999999584E-2)*(cosqf*9.437600198272378E-2+sinqf*9.955366242634227E-1)*1.45384239240002;
	return G2;
}