#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_M33(const std::vector<double>& qs){
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

	double M33 = pow(cosqe*sintheta1*3.333333333333144E-1+cos(qf+9.451666500000044E-2)*costheta1*sinqe*3.333333333333144E-1,2.0)*3.648613099999807E-1+pow(sinqe*sintheta1*3.333333333333144E-1-cos(qf+9.451666500000044E-2)*cosqe*costheta1*3.333333333333144E-1,2.0)*3.648613099999807E-1+(pow(l1*1.380367865060362E1-costheta1*2.162940995833196+l2*costheta1*costheta2*2.927985206164976+l3*costheta1*costheta3*3.973854119136831-l2*sintheta1*sintheta2*5.855970412329953-l3*sintheta1*sintheta3*7.947708238273663,2.0)*4.745645055855675)/(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))+pow(sin(qf+9.451666500000044E-2),2.0)*pow(costheta1,2.0)*1.887401855555595E-1+pow(sinqe*sintheta1-cosqe*cosqf*costheta1*9.955366242634227E-1+cosqe*costheta1*sinqf*9.437600198272378E-2,2.0)*1.482000400000061E-1+pow(cosqe*sintheta1-costheta1*sinqe*sinqf*9.437600198272378E-2+cosqf*costheta1*sinqe*9.955366242634227E-1,2.0)*1.482000400000061E-1+(pow(sintheta1*5.710170987081931E-3-l2*costheta1*sintheta2*1.833053799201778E-2-l2*costheta2*sintheta1*3.642956652745255E-2-l3*costheta1*sintheta3*1.809902853543477E-2-l3*costheta3*sintheta1*3.642956652745255E-2+l2*l3*costheta1*costheta2*sintheta3*1.743114854953092E-1+l2*l3*costheta1*costheta3*sintheta2*1.743114854953092E-1+l2*l3*costheta2*costheta3*sintheta1*1.743114854953092E-1,2.0)*2.514625019249434E1)/((pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0)))+(pow(sintheta1*6.52675530400586E-2-l2*costheta1*sintheta2*2.095190079865006E-1-l2*costheta2*sintheta1*4.163918507756534E-1-l3*costheta1*sintheta3*2.068728427891529E-1-l3*costheta3*sintheta1*4.163918507756534E-1+l2*l3*costheta1*costheta2*sintheta3*1.992389396183398+l2*l3*costheta1*costheta3*sintheta2*1.992389396183398+l2*l3*costheta2*costheta3*sintheta1*1.992389396183398,2.0)*5.536510780289973E1)/((pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0)));
	return M33;
}