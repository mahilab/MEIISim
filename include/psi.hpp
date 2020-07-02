#pragma once
#include <Mahi/Util.hpp>
#include <Eigen/Dense>
using namespace mahi::util;

inline Eigen::VectorXd get_psi(const std::vector<double>& qs){
	Eigen::VectorXd psi = Eigen::VectorXd::Zero(14); 

	double qe = qs[0];
	double qf = qs[1];
	double l1 = qs[2];
	double l2 = qs[3];
	double l3 = qs[4];
	double theta1 = qs[5];
	double theta2 = qs[6];
	double theta3 = qs[7];
	double P_p_x = qs[8];
	double P_p_y = qs[9];
	double P_p_z = qs[10];
	double R_p_x = qs[11];
	double R_p_y = qs[12];
	double R_p_z = qs[13];
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
	double sinR_p_x = sin(R_p_x);
	double cosR_p_x = cos(R_p_x);
	double sinR_p_y = sin(R_p_y);
	double cosR_p_y = cos(R_p_y);
	double sinR_p_z = sin(R_p_z);
	double cosR_p_z = cos(R_p_z);

	psi(0,0) = P_p_x*-1.0+l1*sintheta1-cosR_p_x*sinR_p_y*4.608947781569839E-3-sinR_p_x*sinR_p_y*5.268051420404163E-2+cosR_p_x*cosR_p_y*sinR_p_z*5.268051420404163E-2-cosR_p_y*sinR_p_x*sinR_p_z*4.608947781569839E-3;
	psi(1,0) = P_p_y*-1.0-l1*costheta1*9.955366242634227E-1-cosR_p_x*cosR_p_z*5.268051420404163E-2+cosR_p_z*sinR_p_x*4.608947781569839E-3+1.04065380633557E-1;
	psi(2,0) = P_p_z*-1.0-l1*costheta1*9.437600198272378E-2-cosR_p_x*cosR_p_y*4.608947781569839E-3-cosR_p_y*sinR_p_x*5.268051420404163E-2-cosR_p_x*sinR_p_y*sinR_p_z*5.268051420404163E-2+sinR_p_x*sinR_p_y*sinR_p_z*4.608947781569839E-3+9.480188211044904E-3;
	psi(3,0) = P_p_x*-1.0+l2*sintheta2+cosR_p_x*sinR_p_y*4.792713747590938E-2+sinR_p_x*sinR_p_y*2.234879123846412E-2-cosR_p_x*cosR_p_y*sinR_p_z*2.234879123846412E-2+cosR_p_y*sinR_p_x*sinR_p_z*4.792713747590938E-2;
	psi(4,0) = P_p_y*-1.0+l2*costheta2*4.160362969070661E-1+cosR_p_x*cosR_p_z*2.234879123846412E-2-cosR_p_z*sinR_p_x*4.792713747590938E-2-4.382260649335734E-2;
	psi(5,0) = P_p_z*-1.0+l2*costheta2*9.093480080013023E-1+cosR_p_x*cosR_p_y*4.792713747590938E-2+cosR_p_y*sinR_p_x*2.234879123846412E-2+cosR_p_x*sinR_p_y*sinR_p_z*2.234879123846412E-2-sinR_p_x*sinR_p_y*sinR_p_z*4.792713747590938E-2-9.486335738867524E-2;
	psi(6,0) = P_p_x*-1.0+l3*sintheta3-cosR_p_x*sinR_p_y*4.331818969433954E-2+sinR_p_x*sinR_p_y*3.033172296557396E-2-cosR_p_x*cosR_p_y*sinR_p_z*3.033172296557396E-2-cosR_p_y*sinR_p_x*sinR_p_z*4.331818969433954E-2;
	psi(7,0) = P_p_y*-1.0+l3*costheta3*5.795003273562997E-1+cosR_p_x*cosR_p_z*3.033172296557396E-2+cosR_p_z*sinR_p_x*4.331818969433954E-2-6.024277414019963E-2;
	psi(8,0) = P_p_z*-1.0-l3*costheta3*8.149720060185928E-1-cosR_p_x*cosR_p_y*4.331818969433954E-2+cosR_p_y*sinR_p_x*3.033172296557396E-2+cosR_p_x*sinR_p_y*sinR_p_z*3.033172296557396E-2+sinR_p_x*sinR_p_y*sinR_p_z*4.331818969433954E-2+8.53831691776179E-2;
	psi(9,0) = qe;
	psi(10,0) = qf;
	psi(11,0) = l1;
	psi(12,0) = l2;
	psi(13,0) = l3;
	return psi;
}