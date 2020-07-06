function [qp] = MEII_IK(q_ser)
% parameters
R = 0.1044956; % [m]
r = 0.05288174521; % [m]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]
alpha_5 = 0.094516665; % [rad]
alpha_13 = 5*pi/180; % [rad]

a_ = q_ser(1);
b_ = q_ser(2);
xc = q_ser(3);

c_ = 2*atan((r*cos(a_)*cos(alpha_5)*cos(alpha_13) - (r^2 - 4*a56^2 + r^2*cos(a_)^2*cos(b_)^2 + 2*r^2*cos(a_)*cos(b_))^(1/2) + r*cos(alpha_5)*cos(alpha_13)*cos(b_) + r*cos(a_)*sin(alpha_5)*sin(alpha_13) + r*cos(b_)*sin(alpha_5)*sin(alpha_13) - r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_) + r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_))/(r*cos(a_)*cos(alpha_5)*sin(alpha_13) - 2*a56 - r*cos(a_)*cos(alpha_13)*sin(alpha_5) + r*cos(alpha_5)*cos(b_)*sin(alpha_13) - r*cos(alpha_13)*cos(b_)*sin(alpha_5) + r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_) + r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)));
yc = (3*r*cos(a_)*cos(alpha_13)*cos(c_))/2 - 2*a56*sin(alpha_5) + (r*cos(alpha_13)*cos(b_)*cos(c_))/2 - (3*r*cos(a_)*sin(alpha_13)*sin(c_))/2 - (r*cos(b_)*sin(alpha_13)*sin(c_))/2 - (3*r*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_))/2 - (3*r*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_))/2 - 2*r*cos(a_)*cos(alpha_5)^2*cos(alpha_13)*cos(c_) + 2*r*cos(a_)*cos(alpha_5)^2*sin(alpha_13)*sin(c_) + 2*r*cos(alpha_5)^2*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 2*r*cos(alpha_5)^2*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - 2*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) - 2*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13);
zc = 2*a56*cos(alpha_5) + (r*cos(a_)*cos(alpha_13)*sin(c_))/2 + (r*cos(a_)*cos(c_)*sin(alpha_13))/2 - (r*cos(alpha_13)*cos(b_)*sin(c_))/2 - (r*cos(b_)*cos(c_)*sin(alpha_13))/2 + (r*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_))/2 - (r*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_))/2 + 2*r*cos(alpha_5)^2*cos(alpha_13)*cos(b_)*sin(c_) + 2*r*cos(alpha_5)^2*cos(b_)*cos(c_)*sin(alpha_13) - 2*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(alpha_5) + 2*r*cos(a_)*cos(alpha_5)*sin(alpha_5)*sin(alpha_13)*sin(c_) + 2*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) + 2*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_);

P1 = [0; R*cos(alpha_5) + a56*sin(alpha_5); R*sin(alpha_5) - a56*cos(alpha_5)];
P2 = [0; - R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2);   a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2)];
P3 = [0; - R*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2) - a56*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2);   a56*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2) - R*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2)];

B1 = [xc + r*cos(alpha_13)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)) + r*sin(alpha_13)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)); yc + r*cos(alpha_13)*cos(b_)*cos(c_) - r*cos(b_)*sin(alpha_13)*sin(c_); zc + r*cos(alpha_13)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) + r*sin(alpha_13)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];
B2 = [xc - r*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)) - r*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)); yc - r*cos(b_)*cos(c_)*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2) + r*cos(b_)*sin(c_)*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2); zc - r*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) - r*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];
B3 = [xc - r*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)) - r*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)); yc - r*cos(b_)*cos(c_)*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2) + r*cos(b_)*sin(c_)*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2); zc - r*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) - r*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];

q1 = norm(P1-B1);
q2 = norm(P2-B2);
q3 = norm(P3-B3);

theta1 = atan2(abs(B1(1) - P1(1)),norm(B1(2:3) - P1(2:3)));
theta2 = atan2(abs(B2(1) - P2(1)),norm(B2(2:3) - P2(2:3)));
theta3 = atan2(abs(B3(1) - P3(1)),norm(B3(2:3) - P3(2:3)));

qp = [theta1, theta2, theta3, q1, q2, q3, xc, yc, zc, a_, b_, c_].';

end