#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_M81(const std::vector<double>& qs){
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

	double M81 = cosqf*-1.721869322673653E-4-sinqf*sintheta3*2.421526321482736E-4-(l3*costheta3*sinqe*-3.333333333333144E-1+l3*cos(qf+9.451666500000044E-2)*cosqe*sintheta3*1.666666666666572E-1+l3*sin(qf+9.451666500000044E-2)*cosqe*sintheta3*2.886751345947687E-1)*(cosqe*1.593849999999861E-1+l1*cosqe*sintheta1*3.333333333333144E-1+l2*cosqe*sintheta2*3.333333333333144E-1+l3*cosqe*sintheta3*3.333333333333144E-1+l1*cos(qf+9.451666500000044E-2)*costheta1*sinqe*3.333333333333144E-1-l2*cos(qf+9.451666500000044E-2)*costheta2*sinqe*1.666666666666572E-1-l3*cos(qf+9.451666500000044E-2)*costheta3*sinqe*1.666666666666572E-1+l2*sin(qf+9.451666500000044E-2)*costheta2*sinqe*2.886751345947687E-1-l3*sin(qf+9.451666500000044E-2)*costheta3*sinqe*2.886751345947687E-1)*3.648613099999807E-1-(l3*cosqe*costheta3*3.333333333333144E-1+l3*cos(qf+9.451666500000044E-2)*sinqe*sintheta3*1.666666666666572E-1+l3*sin(qf+9.451666500000044E-2)*sinqe*sintheta3*2.886751345947687E-1)*(sinqe*1.593849999999861E-1+l1*sinqe*sintheta1*3.333333333333144E-1+l2*sinqe*sintheta2*3.333333333333144E-1+l3*sinqe*sintheta3*3.333333333333144E-1-l1*cos(qf+9.451666500000044E-2)*cosqe*costheta1*3.333333333333144E-1+l2*cos(qf+9.451666500000044E-2)*cosqe*costheta2*1.666666666666572E-1+l3*cos(qf+9.451666500000044E-2)*cosqe*costheta3*1.666666666666572E-1-l2*sin(qf+9.451666500000044E-2)*cosqe*costheta2*2.886751345947687E-1+l3*sin(qf+9.451666500000044E-2)*cosqe*costheta3*2.886751345947687E-1)*3.648613099999807E-1+(cosqf*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l1*sintheta1*5.688112926626587E-2-l2*sintheta2*5.691801443319946E-1+l3*sintheta3*5.122990150657643E-1+l1*l2*costheta1*sintheta2*3.44864002803979+l1*l2*costheta2*sintheta1*1.441192008071766-l1*l3*costheta1*sintheta3*3.44864002803979-l1*l3*costheta3*sintheta1*2.007448019968251-l2*l3*costheta2*sintheta3*1.441192008071766+l2*l3*costheta3*sintheta2*2.007448019968251)*1.2E1+sinqf*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l1*sintheta1*-6.24392283801285E-1+l2*sintheta2*2.62935638960073E-1+l3*sintheta3*3.61456644841212E-1+l1*l2*costheta1*sintheta2*3.26928060898581E-1+l1*l2*costheta2*sintheta1*3.150073903239445-l1*l3*costheta1*sintheta3*3.26928060898581E-1+l1*l3*costheta3*sintheta1*2.823145842341091-l2*l3*costheta2*sintheta3*3.150073903239445-l2*l3*costheta3*sintheta2*2.823145842341091)*1.2E1)*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l3*sintheta3*6.261392430703836E-1-l1*l3*costheta1*sintheta3*3.450919662650904-l1*l3*costheta3*sintheta1*6.901839325301808+l2*l3*costheta2*sintheta3*1.463992603082488+l2*l3*costheta3*sintheta2*2.927985206164976)*6.561729075690437E-2-(l3-8.803883999999584E-2)*(costheta3*sinqe*-1.0+cos(qf+9.451666500000044E-2)*cosqe*sintheta3*5.0E-1+sin(qf+9.451666500000044E-2)*cosqe*sintheta3*8.660254037844197E-1)*(cosqe*1.593849999999861E-1-(l3-8.803883999999584E-2)*(cosqe*sintheta3*-1.0+cos(qf+9.451666500000044E-2)*costheta3*sinqe*5.0E-1+sin(qf+9.451666500000044E-2)*costheta3*sinqe*8.660254037844197E-1)*1.0-cos(qf+6.18115440598217E-1)*sinqe*3.833999999999782E-4+sin(qf+6.18115440598217E-1)*sinqe*1.04495600000007E-1)*1.482000400000061E-1-(l3-8.803883999999584E-2)*(cosqe*costheta3+cos(qf+9.451666500000044E-2)*sinqe*sintheta3*5.0E-1+sin(qf+9.451666500000044E-2)*sinqe*sintheta3*8.660254037844197E-1)*(sinqe*1.593849999999861E-1+(l3-8.803883999999584E-2)*(sinqe*sintheta3+cos(qf+9.451666500000044E-2)*cosqe*costheta3*5.0E-1+sin(qf+9.451666500000044E-2)*cosqe*costheta3*8.660254037844197E-1)+cos(qf+6.18115440598217E-1)*cosqe*3.833999999999782E-4-sin(qf+6.18115440598217E-1)*cosqe*1.04495600000007E-1)*1.482000400000061E-1-(cosqf*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l1*costheta1*3.760674912091417E-1+l2*costheta2*1.5372282978351-l3*costheta3*1.869794955751786-2.245306135918668E-3)*1.732050807568839-sinqf*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l1*costheta1*3.966993227389139+l2*costheta2*7.03298146479483E-1+l3*costheta3*1.329550930517144-6.269737997217817E-1)*1.732050807568839)*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l3*costheta3*5.36682734819891E-2-(l3*l3)*3.423911374368345E-1+l1*(l3*l3)*costheta1*1.638304088577797+l2*(l3*l3)*costheta2*1.638304088577797+l1*l3*sintheta1*sintheta3*1.722835144960015E-1+l2*l3*sintheta2*sintheta3*1.701076229408329E-1-l1*l3*costheta1*costheta3*3.423911374368345E-1-l2*l3*costheta2*costheta3*3.423911374368345E-1+l1*l2*l3*costheta1*costheta2*costheta3*1.638304088577797-l1*l2*l3*costheta1*sintheta2*sintheta3*1.638304088577797-l1*l2*l3*costheta2*sintheta1*sintheta3*1.638304088577797)*1.108147996388738E-1-(cosqf*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(l1*costheta1*3.290164220143055E-2+l2*costheta2*3.296596723110724+l3*costheta3*2.670343939073518-6.269737997217817E-1)*1.732050807568839-sinqf*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*(l1*costheta1*3.470669356806866E-1+l2*costheta2*1.508227742306644-l3*costheta3*1.898795511280241+2.245306135918668E-3)*1.732050807568839)*1.0/sqrt(pow(fabs(l1*sintheta1*3.019163492206758E-1-l2*sintheta2*3.139542268885634+l3*sintheta3*2.837625919664788),2.0)*4.0+pow(fabs(l1*costheta1*5.698731594533513E-2+l2*costheta2*5.709873016493475+l3*costheta3*4.625171376159415-1.085950476132211),2.0)+pow(fabs(l1*costheta1*6.011375662261571E-1+l2*costheta2*2.612327079059924-l3*costheta3*3.288810298721273+3.888984305959298E-3),2.0))*1.0/sqrt(pow(fabs(l1*costheta1*6.871033823119433+l2*costheta2*1.21814812257162+l3*costheta3*2.30284976290568-1.085950476132211),2.0)+pow(fabs(l1*sintheta1*-3.450919662650904+l2*sintheta2*1.463992603082488+l3*sintheta3*1.986927059568416),2.0)*4.0+pow(fabs(l1*costheta1*-6.513680018491641E-1-l2*costheta2*2.662557514683158+l3*costheta3*3.238579863098039+3.888984305959298E-3),2.0))*(l3*costheta3*3.757892965479215E-2-(l3*l3)*2.397448554494588E-1+l1*(l3*l3)*costheta1*1.147152872702009+l2*(l3*l3)*costheta2*1.147152872702009+l1*l3*sintheta1*sintheta3*1.206342155593632E-1+l2*l3*sintheta2*sintheta3*1.191106398900956E-1-l1*l3*costheta1*costheta3*2.397448554494588E-1-l2*l3*costheta2*costheta3*2.397448554494588E-1+l1*l2*l3*costheta1*costheta2*costheta3*1.147152872702009-l1*l2*l3*costheta1*sintheta2*sintheta3*1.147152872702009-l1*l2*l3*costheta2*sintheta1*sintheta3*1.147152872702009)*2.439836270297633E-1;
	return M81;
}