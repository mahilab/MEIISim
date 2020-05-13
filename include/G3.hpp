#pragma once
#include <Mahi/Util.hpp>
using namespace mahi::util;

inline double get_G3(const std::vector<double>& qs){
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

	double G3 = cos(qf+9.451666499990097E-2)*costheta2*7.269211961975088E-1-sin(qf+9.451666499990097E-2)*costheta2*1.259064444922842+costheta2*(cosqf*4.160362969123526E-1-sinqf*9.093480079900473E-1)*1.193096483708359-((cos(qf+9.451666499990097E-2)*sintheta2*-6.523679804786298E-3+sin(qf+9.451666499990097E-2)*sintheta2*1.129934487425999E-2+l2*cos(qf+9.451666499990097E-2)*sintheta2*7.410001999960514E-2-l2*sin(qf+9.451666499990097E-2)*sintheta2*1.283449994807597E-1+l2*sintheta2*(cosqf*4.160362969123526E-1-sinqf*9.093480079900473E-1)*1.216204366664897E-1)*(costheta2*sintheta1*sintheta3*-3.080800829229702E-2+l2*sintheta1*sintheta3*1.96548829520907E-1-l1*l2*costheta3*pow(sintheta1,2.0)*6.283017365494743E-1-l2*l3*costheta1*pow(sintheta3,2.0)*6.256454634421971E-1+l1*costheta2*costheta3*pow(sintheta1,2.0)*9.869076803988719E-2+l3*costheta1*costheta2*pow(sintheta3,2.0)*9.785806148101983E-2+l2*pow(costheta2,2.0)*sintheta1*sintheta3*9.827441476045351E-2-(l2*l2)*costheta1*sintheta2*sintheta3*6.283017365494743E-1-(l2*l2)*costheta2*sintheta1*sintheta3*6.269735999958357E-1-(l2*l2)*costheta3*sintheta1*sintheta2*6.256454634421971E-1-l1*pow(sintheta1,2.0)*sintheta2*sintheta3*9.785806148101983E-2-l3*sintheta1*sintheta2*pow(sintheta3,2.0)*9.869076803988719E-2+l1*costheta1*costheta2*sintheta1*sintheta3*9.827441476045351E-2+l2*costheta1*costheta2*sintheta2*sintheta3*9.869076803988719E-2+l2*costheta2*costheta3*sintheta1*sintheta2*9.785806148101983E-2+l3*costheta2*costheta3*sintheta1*sintheta3*9.827441476045351E-2-l1*l2*pow(costheta2,2.0)*costheta3*pow(sintheta1,2.0)*3.141508682747371E-1+l1*(l2*l2)*costheta2*costheta3*pow(sintheta1,2.0)*2.0-l1*l3*costheta2*pow(costheta3,2.0)*pow(sintheta1,2.0)*3.141508682747371E-1-l1*l3*pow(costheta1,2.0)*costheta2*pow(sintheta3,2.0)*3.128227317210985E-1-l2*l3*costheta1*pow(costheta2,2.0)*pow(sintheta3,2.0)*3.128227317210985E-1+(l2*l2)*l3*costheta1*costheta2*pow(sintheta3,2.0)*2.0+l1*l2*costheta3*pow(sintheta1,2.0)*pow(sintheta2,2.0)*6.256454634421971E-1+l1*(l2*l2)*pow(costheta1,2.0)*sintheta2*sintheta3*2.0+l2*l3*costheta1*pow(sintheta2,2.0)*pow(sintheta3,2.0)*6.283017365494743E-1+(l2*l2)*l3*pow(costheta3,2.0)*sintheta1*sintheta2*2.0-l1*l2*costheta1*sintheta1*sintheta3*3.141508682747371E-1-l2*l3*costheta3*sintheta1*sintheta3*3.128227317210985E-1-l1*l2*costheta1*pow(costheta2,2.0)*sintheta1*sintheta3*3.134867999979178E-1+l1*(l2*l2)*costheta1*costheta2*sintheta1*sintheta3+l1*(l2*l2)*costheta1*costheta3*sintheta1*sintheta2*4.0-l1*l2*pow(costheta1,2.0)*costheta2*sintheta2*sintheta3*3.141508682747371E-1-l2*l3*costheta2*pow(costheta3,2.0)*sintheta1*sintheta2*3.128227317210985E-1-l2*l3*pow(costheta2,2.0)*costheta3*sintheta1*sintheta3*3.134867999979178E-1+(l2*l2)*l3*costheta1*costheta3*sintheta2*sintheta3*4.0+(l2*l2)*l3*costheta2*costheta3*sintheta1*sintheta3+l1*l2*costheta2*pow(sintheta1,2.0)*sintheta2*sintheta3*3.128227317210985E-1+l1*l3*costheta1*sintheta1*sintheta2*pow(sintheta3,2.0)*9.397963317169342E-1+l1*l3*costheta3*pow(sintheta1,2.0)*sintheta2*sintheta3*9.411244682851247E-1+l2*l3*costheta2*sintheta1*sintheta2*pow(sintheta3,2.0)*3.141508682747371E-1+l1*l2*l3*pow(costheta1,2.0)*pow(costheta2,2.0)*pow(sintheta3,2.0)+l1*l2*l3*pow(costheta2,2.0)*pow(costheta3,2.0)*pow(sintheta1,2.0)-l1*l2*l3*pow(costheta1,2.0)*pow(sintheta2,2.0)*pow(sintheta3,2.0)*2.0-l1*l2*l3*pow(costheta3,2.0)*pow(sintheta1,2.0)*pow(sintheta2,2.0)*2.0-l1*l2*costheta1*costheta2*costheta3*sintheta1*sintheta2*9.411244682851247E-1-l1*l3*costheta1*costheta2*costheta3*sintheta1*sintheta3*3.134867999979178E-1-l2*l3*costheta1*costheta2*costheta3*sintheta2*sintheta3*9.397963317169342E-1+l1*l2*l3*costheta1*costheta2*pow(costheta3,2.0)*sintheta1*sintheta2*3.0+l1*l2*l3*costheta1*pow(costheta2,2.0)*costheta3*sintheta1*sintheta3+l1*l2*l3*pow(costheta1,2.0)*costheta2*costheta3*sintheta2*sintheta3*3.0-l1*l2*l3*costheta1*costheta2*sintheta1*sintheta2*pow(sintheta3,2.0)*3.0-l1*l2*l3*costheta1*costheta3*sintheta1*pow(sintheta2,2.0)*sintheta3*8.0-l1*l2*l3*costheta2*costheta3*pow(sintheta1,2.0)*sintheta2*sintheta3*3.0)*9.810000000055879)/((l2*l2)*costheta1*pow(sintheta2,2.0)*sintheta3*9.869076803988719E-2+(l2*l2)*costheta3*sintheta1*pow(sintheta2,2.0)*9.785806148101983E-2-l2*sintheta1*sintheta2*sintheta3*3.080800829229702E-2+l1*l2*costheta2*pow(sintheta1,2.0)*sintheta3*9.785806148101983E-2+l1*l2*costheta3*pow(sintheta1,2.0)*sintheta2*9.869076803988719E-2+l2*l3*costheta1*sintheta2*pow(sintheta3,2.0)*9.785806148101983E-2+l2*l3*costheta2*sintheta1*pow(sintheta3,2.0)*9.869076803988719E-2+(l2*l2)*costheta2*sintheta1*sintheta2*sintheta3*9.827441476045351E-2-l1*(l2*l2)*pow(costheta1,2.0)*pow(sintheta2,2.0)*sintheta3*3.141508682747371E-1-l1*(l2*l2)*pow(costheta2,2.0)*pow(sintheta1,2.0)*sintheta3*3.128227317210985E-1-(l2*l2)*l3*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.141508682747371E-1-(l2*l2)*l3*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.128227317210985E-1-l1*l2*l3*pow(costheta1,2.0)*sintheta2*pow(sintheta3,2.0)*3.128227317210985E-1-l1*l2*l3*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.141508682747371E-1-l1*(l2*l2)*costheta1*costheta3*sintheta1*pow(sintheta2,2.0)*9.411244682851247E-1-l1*(l2*l2)*costheta2*costheta3*pow(sintheta1,2.0)*sintheta2*9.397963317169342E-1-(l2*l2)*l3*costheta1*costheta2*sintheta2*pow(sintheta3,2.0)*9.411244682851247E-1-(l2*l2)*l3*costheta1*costheta3*pow(sintheta2,2.0)*sintheta3*9.397963317169342E-1+l1*l2*costheta1*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+l2*l3*costheta3*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+l1*(l2*l2)*l3*costheta1*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.0+l1*(l2*l2)*l3*costheta1*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.0+l1*(l2*l2)*l3*costheta2*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.0+l1*(l2*l2)*l3*pow(costheta1,2.0)*costheta2*sintheta2*pow(sintheta3,2.0)*3.0+l1*(l2*l2)*l3*pow(costheta1,2.0)*costheta3*pow(sintheta2,2.0)*sintheta3*3.0+l1*(l2*l2)*l3*pow(costheta2,2.0)*costheta3*pow(sintheta1,2.0)*sintheta3*3.0-l1*l2*l3*costheta1*costheta2*sintheta1*pow(sintheta3,2.0)*9.397963317169342E-1-l1*l2*l3*costheta2*costheta3*pow(sintheta1,2.0)*sintheta3*9.411244682851247E-1-l1*(l2*l2)*costheta1*costheta2*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-(l2*l2)*l3*costheta2*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-l1*l2*l3*costheta1*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1+l1*(l2*l2)*l3*costheta1*costheta2*costheta3*sintheta1*sintheta2*sintheta3*9.0)-((cos(qf+9.451666499990097E-2)*sintheta3*-6.523679804786298E-3-sin(qf+9.451666499990097E-2)*sintheta3*1.129934487425999E-2+l3*sintheta3*(cosqf*5.795003273524344E-1+sinqf*8.149720060173422E-1)*1.216204366664897E-1+l3*cos(qf+9.451666499990097E-2)*sintheta3*7.410001999960514E-2+l3*sin(qf+9.451666499990097E-2)*sintheta3*1.283449994807597E-1)*(l1*9.827353278706141E-2-l1*pow(costheta1,2.0)*9.827353278706141E-2+(l3*l3)*costheta1*6.256454634421971E-1-l1*l2*costheta2*6.256454634421971E-1-l1*l3*costheta3*6.269735999958357E-1-(l3*l3)*costheta1*pow(costheta3,2.0)*6.256454634421971E-1-l2*sintheta1*sintheta2*4.154713065602778E-4-l3*sintheta1*sintheta3*9.785806148101983E-2-l1*(l3*l3)*pow(costheta1,2.0)*2.0+l1*(l3*l3)*pow(costheta3,2.0)+l1*(l3*l3)*pow(costheta1,2.0)*pow(costheta3,2.0)+l1*l2*pow(costheta1,2.0)*costheta2*6.256454634421971E-1+l1*l3*pow(costheta1,2.0)*costheta3*6.269735999958357E-1-l2*(l3*l3)*costheta1*costheta2*4.0+(l3*l3)*costheta3*sintheta1*sintheta3*3.128227317210985E-1+l1*l2*l3*costheta2*costheta3*2.0+l2*(l3*l3)*costheta1*costheta2*pow(costheta3,2.0)*4.0-l2*(l3*l3)*pow(costheta3,2.0)*sintheta1*sintheta2*1.0-l1*l2*costheta1*sintheta1*sintheta2*3.128227317210985E-1-l1*l3*costheta1*sintheta1*sintheta3*3.154790048356517E-1+l2*l3*costheta1*sintheta2*sintheta3*2.656273118475383E-3+l2*l3*costheta2*sintheta1*sintheta3*6.256454634421971E-1+l2*l3*costheta3*sintheta1*sintheta2*3.141508682747371E-1+l1*(l3*l3)*costheta1*costheta3*sintheta1*sintheta3-l2*(l3*l3)*costheta1*costheta3*sintheta2*sintheta3*2.0-l2*(l3*l3)*costheta2*costheta3*sintheta1*sintheta3*2.0-l1*l2*l3*pow(costheta1,2.0)*costheta2*costheta3*2.0+l1*l2*l3*pow(costheta1,2.0)*sintheta2*sintheta3*2.0+l1*l2*l3*costheta1*costheta2*sintheta1*sintheta3*4.0+l1*l2*l3*costheta1*costheta3*sintheta1*sintheta2)*9.810000000055879)/((l3*l3)*costheta1*sintheta2*pow(sintheta3,2.0)*9.785806148101983E-2+(l3*l3)*costheta2*sintheta1*pow(sintheta3,2.0)*9.869076803988719E-2-l3*sintheta1*sintheta2*sintheta3*3.080800829229702E-2+l1*l3*costheta2*pow(sintheta1,2.0)*sintheta3*9.785806148101983E-2+l1*l3*costheta3*pow(sintheta1,2.0)*sintheta2*9.869076803988719E-2+l2*l3*costheta1*pow(sintheta2,2.0)*sintheta3*9.869076803988719E-2+l2*l3*costheta3*sintheta1*pow(sintheta2,2.0)*9.785806148101983E-2+(l3*l3)*costheta3*sintheta1*sintheta2*sintheta3*9.827441476045351E-2-l1*(l3*l3)*pow(costheta1,2.0)*sintheta2*pow(sintheta3,2.0)*3.128227317210985E-1-l1*(l3*l3)*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.141508682747371E-1-l2*(l3*l3)*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.141508682747371E-1-l2*(l3*l3)*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.128227317210985E-1-l1*l2*l3*pow(costheta1,2.0)*pow(sintheta2,2.0)*sintheta3*3.141508682747371E-1-l1*l2*l3*pow(costheta2,2.0)*pow(sintheta1,2.0)*sintheta3*3.128227317210985E-1-l1*(l3*l3)*costheta1*costheta2*sintheta1*pow(sintheta3,2.0)*9.397963317169342E-1-l1*(l3*l3)*costheta2*costheta3*pow(sintheta1,2.0)*sintheta3*9.411244682851247E-1-l2*(l3*l3)*costheta1*costheta2*sintheta2*pow(sintheta3,2.0)*9.411244682851247E-1-l2*(l3*l3)*costheta1*costheta3*pow(sintheta2,2.0)*sintheta3*9.397963317169342E-1+l1*l3*costheta1*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+l2*l3*costheta2*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+l1*l2*(l3*l3)*costheta1*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.0+l1*l2*(l3*l3)*costheta1*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.0+l1*l2*(l3*l3)*costheta2*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.0+l1*l2*(l3*l3)*pow(costheta1,2.0)*costheta2*sintheta2*pow(sintheta3,2.0)*3.0+l1*l2*(l3*l3)*pow(costheta1,2.0)*costheta3*pow(sintheta2,2.0)*sintheta3*3.0+l1*l2*(l3*l3)*pow(costheta2,2.0)*costheta3*pow(sintheta1,2.0)*sintheta3*3.0-l1*l2*l3*costheta1*costheta3*sintheta1*pow(sintheta2,2.0)*9.411244682851247E-1-l1*l2*l3*costheta2*costheta3*pow(sintheta1,2.0)*sintheta2*9.397963317169342E-1-l1*(l3*l3)*costheta1*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-l2*(l3*l3)*costheta2*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-l1*l2*l3*costheta1*costheta2*sintheta1*sintheta2*sintheta3*3.134867999979178E-1+l1*l2*(l3*l3)*costheta1*costheta2*costheta3*sintheta1*sintheta2*sintheta3*9.0)+((cos(qf+9.451666499990097E-2)*sintheta1*-1.30473596095726E-2+l1*cos(qf+9.451666499990097E-2)*sintheta1*1.482000399992103E-1+l1*sintheta1*(cosqf*9.95536624264787E-1-sinqf*9.437600198180007E-2)*1.216204366664897E-1)*(l3*9.827353278706141E-2+(l1*l1)*costheta3*6.283017365494743E-1-l3*pow(costheta3,2.0)*9.827353278706141E-2-l1*l3*costheta1*6.269735999958357E-1-l2*l3*costheta2*6.283017365494743E-1-(l1*l1)*pow(costheta1,2.0)*costheta3*6.283017365494743E-1-l1*sintheta1*sintheta3*9.869076803988719E-2+l2*sintheta2*sintheta3*4.17235253273418E-4+(l1*l1)*l3*pow(costheta1,2.0)-(l1*l1)*l3*pow(costheta3,2.0)*2.0+(l1*l1)*l3*pow(costheta1,2.0)*pow(costheta3,2.0)+l1*l3*costheta1*pow(costheta3,2.0)*6.269735999958357E-1-(l1*l1)*l2*costheta2*costheta3*4.0+l2*l3*costheta2*pow(costheta3,2.0)*6.283017365494743E-1+(l1*l1)*costheta1*sintheta1*sintheta3*3.141508682747371E-1+l1*l2*l3*costheta1*costheta2*2.0+(l1*l1)*l2*pow(costheta1,2.0)*costheta2*costheta3*4.0-(l1*l1)*l2*pow(costheta1,2.0)*sintheta2*sintheta3*1.0+l1*l2*costheta1*sintheta2*sintheta3*3.128227317210985E-1+l1*l2*costheta2*sintheta1*sintheta3*6.283017365494743E-1-l1*l2*costheta3*sintheta1*sintheta2*2.656273118475383E-3-l1*l3*costheta3*sintheta1*sintheta3*3.11494595160184E-1-l2*l3*costheta3*sintheta2*sintheta3*3.141508682747371E-1-(l1*l1)*l2*costheta1*costheta2*sintheta1*sintheta3*2.0-(l1*l1)*l2*costheta1*costheta3*sintheta1*sintheta2*2.0+(l1*l1)*l3*costheta1*costheta3*sintheta1*sintheta3-l1*l2*l3*costheta1*costheta2*pow(costheta3,2.0)*2.0+l1*l2*l3*pow(costheta3,2.0)*sintheta1*sintheta2*2.0+l1*l2*l3*costheta1*costheta3*sintheta2*sintheta3+l1*l2*l3*costheta2*costheta3*sintheta1*sintheta3*4.0)*9.810000000055879)/((l1*l1)*costheta2*pow(sintheta1,2.0)*sintheta3*9.785806148101983E-2+(l1*l1)*costheta3*pow(sintheta1,2.0)*sintheta2*9.869076803988719E-2-l1*sintheta1*sintheta2*sintheta3*3.080800829229702E-2+l1*l2*costheta1*pow(sintheta2,2.0)*sintheta3*9.869076803988719E-2+l1*l2*costheta3*sintheta1*pow(sintheta2,2.0)*9.785806148101983E-2+l1*l3*costheta1*sintheta2*pow(sintheta3,2.0)*9.785806148101983E-2+l1*l3*costheta2*sintheta1*pow(sintheta3,2.0)*9.869076803988719E-2+(l1*l1)*costheta1*sintheta1*sintheta2*sintheta3*9.827441476045351E-2-(l1*l1)*l2*pow(costheta1,2.0)*pow(sintheta2,2.0)*sintheta3*3.141508682747371E-1-(l1*l1)*l2*pow(costheta2,2.0)*pow(sintheta1,2.0)*sintheta3*3.128227317210985E-1-(l1*l1)*l3*pow(costheta1,2.0)*sintheta2*pow(sintheta3,2.0)*3.128227317210985E-1-(l1*l1)*l3*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.141508682747371E-1-l1*l2*l3*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.141508682747371E-1-l1*l2*l3*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.128227317210985E-1-(l1*l1)*l2*costheta1*costheta3*sintheta1*pow(sintheta2,2.0)*9.411244682851247E-1-(l1*l1)*l2*costheta2*costheta3*pow(sintheta1,2.0)*sintheta2*9.397963317169342E-1-(l1*l1)*l3*costheta1*costheta2*sintheta1*pow(sintheta3,2.0)*9.397963317169342E-1-(l1*l1)*l3*costheta2*costheta3*pow(sintheta1,2.0)*sintheta3*9.411244682851247E-1+l1*l2*costheta2*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+l1*l3*costheta3*sintheta1*sintheta2*sintheta3*9.827441476045351E-2+(l1*l1)*l2*l3*costheta1*pow(costheta2,2.0)*sintheta1*pow(sintheta3,2.0)*3.0+(l1*l1)*l2*l3*costheta1*pow(costheta3,2.0)*sintheta1*pow(sintheta2,2.0)*3.0+(l1*l1)*l2*l3*costheta2*pow(costheta3,2.0)*pow(sintheta1,2.0)*sintheta2*3.0+(l1*l1)*l2*l3*pow(costheta1,2.0)*costheta2*sintheta2*pow(sintheta3,2.0)*3.0+(l1*l1)*l2*l3*pow(costheta1,2.0)*costheta3*pow(sintheta2,2.0)*sintheta3*3.0+(l1*l1)*l2*l3*pow(costheta2,2.0)*costheta3*pow(sintheta1,2.0)*sintheta3*3.0-l1*l2*l3*costheta1*costheta2*sintheta2*pow(sintheta3,2.0)*9.411244682851247E-1-l1*l2*l3*costheta1*costheta3*pow(sintheta2,2.0)*sintheta3*9.397963317169342E-1-(l1*l1)*l2*costheta1*costheta2*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-(l1*l1)*l3*costheta1*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1-l1*l2*l3*costheta2*costheta3*sintheta1*sintheta2*sintheta3*3.134867999979178E-1+(l1*l1)*l2*l3*costheta1*costheta2*costheta3*sintheta1*sintheta2*sintheta3*9.0);
	return G3;
}