#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_M23(const std::vector<double>& qs){
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

	double M23 = (sinqe*sintheta1*3.333333333333144E-1-cos(qf+9.451666500000044E-2)*cosqe*costheta1*3.333333333333144E-1)*(l2*cos(qf+9.451666500000044E-2)*cosqe*costheta2*2.886751345947687E-1-l3*cos(qf+9.451666500000044E-2)*cosqe*costheta3*2.886751345947687E-1-l1*sin(qf+9.451666500000044E-2)*cosqe*costheta1*3.333333333333144E-1+l2*sin(qf+9.451666500000044E-2)*cosqe*costheta2*1.666666666666572E-1+l3*sin(qf+9.451666500000044E-2)*cosqe*costheta3*1.666666666666572E-1)*-3.648613099999807E-1+(cosqe*sintheta1*3.333333333333144E-1+cos(qf+9.451666500000044E-2)*costheta1*sinqe*3.333333333333144E-1)*(l2*cos(qf+9.451666500000044E-2)*costheta2*sinqe*2.886751345947687E-1-l3*cos(qf+9.451666500000044E-2)*costheta3*sinqe*2.886751345947687E-1-l1*sin(qf+9.451666500000044E-2)*costheta1*sinqe*3.333333333333144E-1+l2*sin(qf+9.451666500000044E-2)*costheta2*sinqe*1.666666666666572E-1+l3*sin(qf+9.451666500000044E-2)*costheta3*sinqe*1.666666666666572E-1)*3.648613099999807E-1-sin(qf+9.451666500000044E-2)*costheta1*(l1*cos(qf+9.451666500000044E-2)*costheta1*-3.333333333333144E-1+l2*cos(qf+9.451666500000044E-2)*costheta2*1.666666666666572E-1+l3*cos(qf+9.451666500000044E-2)*costheta3*1.666666666666572E-1-l2*sin(qf+9.451666500000044E-2)*costheta2*2.886751345947687E-1+l3*sin(qf+9.451666500000044E-2)*costheta3*2.886751345947687E-1)*1.216204366666602E-1+cosqe*(sinqe*sintheta1-cosqe*cosqf*costheta1*9.955366242634227E-1+cosqe*costheta1*sinqf*9.437600198272378E-2)*(cos(qf+9.451666500000044E-2)*3.833999999999782E-4-sin(qf+9.451666500000044E-2)*1.04495600000007E-1-sin(qf+9.451666500000044E-2)*costheta1*8.803883999999584E-2+l1*sin(qf+9.451666500000044E-2)*costheta1)*1.482000400000061E-1-sin(qf+9.451666500000044E-2)*costheta1*(cos(qf+9.451666500000044E-2)*1.04495600000007E-1+sin(qf+9.451666500000044E-2)*3.833999999999782E-4-cos(qf+9.451666500000044E-2)*costheta1*(l1-8.803883999999584E-2)*1.0)*1.482000400000061E-1-sinqe*(cosqe*sintheta1-costheta1*sinqe*sinqf*9.437600198272378E-2+cosqf*costheta1*sinqe*9.955366242634227E-1)*(cos(qf+9.451666500000044E-2)*3.833999999999782E-4-sin(qf+9.451666500000044E-2)*1.04495600000007E-1-sin(qf+9.451666500000044E-2)*costheta1*8.803883999999584E-2+l1*sin(qf+9.451666500000044E-2)*costheta1)*1.482000400000061E-1-(1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l1*sintheta1*1.743114854953092E-1-l2*sintheta2*1.812615574073106+l3*sintheta3*1.638304088577797)*(sintheta1*6.52675530400586E-2-l2*costheta1*sintheta2*2.095190079865006E-1-l2*costheta2*sintheta1*4.163918507756534E-1-l3*costheta1*sintheta3*2.068728427891529E-1-l3*costheta3*sintheta1*4.163918507756534E-1+l2*l3*costheta1*costheta2*sintheta3*1.992389396183398+l2*l3*costheta1*costheta3*sintheta2*1.992389396183398+l2*l3*costheta2*costheta3*sintheta1*1.992389396183398)*8.45184076460896E-1)/(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))-(1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(l1*1.380367865060362E1-costheta1*2.162940995833196+l2*costheta1*costheta2*2.927985206164976+l3*costheta1*costheta3*3.973854119136831-l2*sintheta1*sintheta2*5.855970412329953-l3*sintheta1*sintheta3*7.947708238273663)*(l1*costheta1*-2.08991200000014E-1-l2*costheta2*2.08991200000014E-1-l3*costheta3*2.08991200000014E-1+l1*l2*costheta1*costheta2+l1*l3*costheta1*costheta3+l2*l3*costheta2*costheta3+3.275843224476205E-2)*2.362222467248557)/(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))-(1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(l1*sintheta1*-1.992389396183398+l2*sintheta2*8.452365234813897E-1+l3*sintheta3*1.147152872702009)*(sintheta1*5.710170987081931E-3-l2*costheta1*sintheta2*1.833053799201778E-2-l2*costheta2*sintheta1*3.642956652745255E-2-l3*costheta1*sintheta3*1.809902853543477E-2-l3*costheta3*sintheta1*3.642956652745255E-2+l2*l3*costheta1*costheta2*sintheta3*1.743114854953092E-1+l2*l3*costheta1*costheta3*sintheta2*1.743114854953092E-1+l2*l3*costheta2*costheta3*sintheta1*1.743114854953092E-1)*3.838737264101724E-1)/(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0));
	return M23;
}