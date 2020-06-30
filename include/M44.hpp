#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_M44(const std::vector<double>& qs){
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

	double M44 = pow(cos(qf+9.451666500000044E-2)*costheta2*2.886751345947687E-1+sin(qf+9.451666500000044E-2)*costheta2*1.666666666666572E-1,2.0)*3.648613099999807E-1+pow(sinqe*sintheta2*3.333333333333144E-1+cos(qf+9.451666500000044E-2)*cosqe*costheta2*1.666666666666572E-1-sin(qf+9.451666500000044E-2)*cosqe*costheta2*2.886751345947687E-1,2.0)*3.648613099999807E-1+pow(cosqe*sintheta2*3.333333333333144E-1-cos(qf+9.451666500000044E-2)*costheta2*sinqe*1.666666666666572E-1+sin(qf+9.451666500000044E-2)*costheta2*sinqe*2.886751345947687E-1,2.0)*3.648613099999807E-1+pow(sin(qf+1.141714216196306),2.0)*pow(costheta2,2.0)*1.482000400000061E-1+(pow(l2*2.927985206164976-costheta2*4.553312548463282E-1+l1*costheta1*costheta2*3.450919662650904-l3*costheta2*costheta3*1.986927059568416-l1*sintheta1*sintheta2*6.901839325301808+l3*sintheta2*sintheta3*3.973854119136831,2.0)*4.745645055855675)/(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))+pow(sinqe*sintheta2+cos(qf+9.451666500000044E-2)*cosqe*costheta2*5.0E-1-sin(qf+9.451666500000044E-2)*cosqe*costheta2*8.660254037844197E-1,2.0)*1.482000400000061E-1+pow(cosqe*sintheta2-cos(qf+9.451666500000044E-2)*costheta2*sinqe*5.0E-1+sin(qf+9.451666500000044E-2)*costheta2*sinqe*8.660254037844197E-1,2.0)*1.482000400000061E-1+(pow(sintheta2*2.76886233852629E-2-l1*costheta1*sintheta2*1.766469953261947E-1-l1*costheta2*sintheta1*8.776220289908565E-2-l3*costheta2*sintheta3*8.888479242712322E-2-l3*costheta3*sintheta2*1.766469953261947E-1+l1*l3*costheta1*costheta2*sintheta3*8.452365234813897E-1+l1*l3*costheta1*costheta3*sintheta2*8.452365234813897E-1+l1*l3*costheta2*costheta3*sintheta1*8.452365234813897E-1,2.0)*5.536510780289973E1)/((pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0)))+(pow(sintheta2*5.937844446907548E-2-l1*costheta1*sintheta2*3.788207039642657E-1-l1*costheta2*sintheta1*1.882066514762357E-1-l3*costheta2*sintheta3*1.9061405248803E-1-l3*costheta3*sintheta2*3.788207039642657E-1+l1*l3*costheta1*costheta2*sintheta3*1.812615574073106+l1*l3*costheta1*costheta3*sintheta2*1.812615574073106+l1*l3*costheta2*costheta3*sintheta1*1.812615574073106,2.0)*2.514625019249434E1)/((pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0)));
	return M44;
}