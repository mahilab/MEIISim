function [qp] = MEII_IK_3(q_ser)
% parameters
R = 0.1044956; % [m]
r = 0.05288174521; % [m]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]
alpha_5 = 0.094516665; % [rad]
alpha_13 = 5*pi/180; % [rad]

a_ = q_ser(1);cd 
b_ = q_ser(2);
xc = q_ser(3);

c_ = 2*atan((r*cos(a_)*cos(alpha_5)*cos(alpha_13) - (r^2 - 4*a56^2 + r^2*cos(a_)^2*cos(b_)^2 + 2*r^2*cos(a_)*cos(b_))^(1/2) + r*cos(alpha_5)*cos(alpha_13)*cos(b_) + r*cos(a_)*sin(alpha_5)*sin(alpha_13) + r*cos(b_)*sin(alpha_5)*sin(alpha_13) - r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_) + r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_))/(r*cos(a_)*cos(alpha_5)*sin(alpha_13) - 2*a56 - r*cos(a_)*cos(alpha_13)*sin(alpha_5) + r*cos(alpha_5)*cos(b_)*sin(alpha_13) - r*cos(alpha_13)*cos(b_)*sin(alpha_5) + r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_) + r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)));
yc = (3*r*cos(a_)*cos(alpha_13)*cos(c_))/2 - 2*a56*sin(alpha_5) + (r*cos(alpha_13)*cos(b_)*cos(c_))/2 - (3*r*cos(a_)*sin(alpha_13)*sin(c_))/2 - (r*cos(b_)*sin(alpha_13)*sin(c_))/2 - (3*r*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_))/2 - (3*r*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_))/2 - 2*r*cos(a_)*cos(alpha_5)^2*cos(alpha_13)*cos(c_) + 2*r*cos(a_)*cos(alpha_5)^2*sin(alpha_13)*sin(c_) + 2*r*cos(alpha_5)^2*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 2*r*cos(alpha_5)^2*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - 2*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) - 2*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13);
zc = 2*a56*cos(alpha_5) + (r*cos(a_)*cos(alpha_13)*sin(c_))/2 + (r*cos(a_)*cos(c_)*sin(alpha_13))/2 - (r*cos(alpha_13)*cos(b_)*sin(c_))/2 - (r*cos(b_)*cos(c_)*sin(alpha_13))/2 + (r*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_))/2 - (r*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_))/2 + 2*r*cos(alpha_5)^2*cos(alpha_13)*cos(b_)*sin(c_) + 2*r*cos(alpha_5)^2*cos(b_)*cos(c_)*sin(alpha_13) - 2*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(alpha_5) + 2*r*cos(a_)*cos(alpha_5)*sin(alpha_5)*sin(alpha_13)*sin(c_) + 2*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) + 2*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_);

P{1} = [0; R*cos(alpha_5) + a56*sin(alpha_5); R*sin(alpha_5) - a56*cos(alpha_5)];
P{2} = [0; - R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2);   a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2)];
P{3} = [0; - R*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2) - a56*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2);   a56*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2) - R*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2)];

B{1} = [xc + r*cos(alpha_13)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)) + r*sin(alpha_13)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)); yc + r*cos(alpha_13)*cos(b_)*cos(c_) - r*cos(b_)*sin(alpha_13)*sin(c_); zc + r*cos(alpha_13)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) + r*sin(alpha_13)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];
B{2} = [xc - r*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)) - r*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)); yc - r*cos(b_)*cos(c_)*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2) + r*cos(b_)*sin(c_)*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2); zc - r*(cos(alpha_13)/2 - (3^(1/2)*sin(alpha_13))/2)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) - r*(sin(alpha_13)/2 + (3^(1/2)*cos(alpha_13))/2)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];
B{3} = [xc - r*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2)*(cos(c_)*sin(a_) + cos(a_)*sin(b_)*sin(c_)) - r*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2)*(sin(a_)*sin(c_) - cos(a_)*cos(c_)*sin(b_)); yc - r*cos(b_)*cos(c_)*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2) + r*cos(b_)*sin(c_)*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2); zc - r*(cos(alpha_13)/2 + (3^(1/2)*sin(alpha_13))/2)*(cos(a_)*sin(c_) + cos(c_)*sin(a_)*sin(b_)) - r*(sin(alpha_13)/2 - (3^(1/2)*cos(alpha_13))/2)*(cos(a_)*cos(c_) - sin(a_)*sin(b_)*sin(c_))];

q1 = sqrt((P{1}(1)-B{1}(1))^2 + (P{1}(2)-B{1}(2))^2 + (P{1}(3)-B{1}(3))^2);
q2 = sqrt((P{2}(1)-B{2}(1))^2 + (P{2}(2)-B{2}(2))^2 + (P{2}(3)-B{2}(3))^2);
q3 = sqrt((P{3}(1)-B{3}(1))^2 + (P{3}(2)-B{3}(2))^2 + (P{3}(3)-B{3}(3))^2);

theta1 = atan2(sqrt((B{1}(1) - P{1}(1))^2),sqrt((B{1}(2) - P{1}(2))^2 + (B{1}(3) - P{1}(3))^2));
theta2 = atan2(sqrt((B{2}(1) - P{2}(1))^2),sqrt((B{2}(2) - P{2}(2))^2 + (B{2}(3) - P{2}(3))^2));
theta3 = atan2(sqrt((B{3}(1) - P{3}(1))^2),sqrt((B{3}(2) - P{3}(2))^2 + (B{3}(3) - P{3}(3))^2));

qp = [theta1, theta2, theta3, q1, q2, q3, xc, yc, zc, a_, b_, c_];

end