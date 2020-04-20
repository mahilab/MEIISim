function q_out = MEII_FK_2(q_in)
% Forward kinematics for the MEII

l1 = q_in(1);
l2 = q_in(2);
l3 = q_in(3);

R = 0.1044956; % [m]
r = 0.05288174521; % [m]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]
alpha_5 = 0.094516665; % [rad]
alpha_13 = 5*pi/180; % [rad]

thetas = [pi/4,pi/4,pi/4].';
old_thetas = [pi pi pi].';
iter = 0;

while abs(rms(old_thetas-thetas)) > 1e-12 && iter < 20
    old_thetas = thetas;
    t1 = thetas(1); t2 = thetas(2); t3 = thetas(3);
    jac = [ (l1*((3*l2*sin(t1 + t2))/2 - 3*R*cos(t1) + (l2*sin(t1 - t2))/2 - 3^(1/2)*a56*cos(t1)))/(2*(3*R^2 + 3*a56^2 + l1^2 + l2^2 - 3*R*l1*sin(t1) - 3*R*l2*sin(t2) - (l1*l2*cos(t1 - t2))/2 - (3*l1*l2*cos(t1 + t2))/2 - 3^(1/2)*a56*l1*sin(t1) + 3^(1/2)*a56*l2*sin(t2))^(1/2)), (l2*((3*l1*sin(t1 + t2))/2 - 3*R*cos(t2) - (l1*sin(t1 - t2))/2 + 3^(1/2)*a56*cos(t2)))/(2*(3*R^2 + 3*a56^2 + l1^2 + l2^2 - 3*R*l1*sin(t1) - 3*R*l2*sin(t2) - (l1*l2*cos(t1 - t2))/2 - (3*l1*l2*cos(t1 + t2))/2 - 3^(1/2)*a56*l1*sin(t1) + 3^(1/2)*a56*l2*sin(t2))^(1/2)),0; 0, (l2*((3*l3*sin(t2 + t3))/2 - 3*R*cos(t2) + (l3*sin(t2 - t3))/2 - 3^(1/2)*a56*cos(t2)))/(2*(3*R^2 + 3*a56^2 + l2^2 + l3^2 - 3*R*l2*sin(t2) - 3*R*l3*sin(t3) - (l2*l3*cos(t2 - t3))/2 - (3*l2*l3*cos(t2 + t3))/2 - 3^(1/2)*a56*l2*sin(t2) + 3^(1/2)*a56*l3*sin(t3))^(1/2)),  (l3*((3*l2*sin(t2 + t3))/2 - 3*R*cos(t3) - (l2*sin(t2 - t3))/2 + 3^(1/2)*a56*cos(t3)))/(2*(3*R^2 + 3*a56^2 + l2^2 + l3^2 - 3*R*l2*sin(t2) - 3*R*l3*sin(t3) - (l2*l3*cos(t2 - t3))/2 - (3*l2*l3*cos(t2 + t3))/2 - 3^(1/2)*a56*l2*sin(t2) + 3^(1/2)*a56*l3*sin(t3))^(1/2)); (l1*((3*l3*sin(t1 + t3))/2 - 3*R*cos(t1) + (l3*sin(t1 - t3))/2 + 3^(1/2)*a56*cos(t1)))/(2*(3*R^2 + 3*a56^2 + l1^2 + l3^2 - 3*R*l1*sin(t1) - 3*R*l3*sin(t3) - (l1*l3*cos(t1 - t3))/2 - (3*l1*l3*cos(t1 + t3))/2 + 3^(1/2)*a56*l1*sin(t1) - 3^(1/2)*a56*l3*sin(t3))^(1/2)),0, -(l3*(3*R*cos(t3) - (3*l1*sin(t1 + t3))/2 + (l1*sin(t1 - t3))/2 + 3^(1/2)*a56*cos(t3)))/(2*(3*R^2 + 3*a56^2 + l1^2 + l3^2 - 3*R*l1*sin(t1) - 3*R*l3*sin(t3) - (l1*l3*cos(t1 - t3))/2 - (3*l1*l3*cos(t1 + t3))/2 + 3^(1/2)*a56*l1*sin(t1) - 3^(1/2)*a56*l3*sin(t3))^(1/2))];
    eq =  [(3*R^2 + 3*a56^2 + l1^2 + l2^2 - 3*R*l1*sin(t1) - 3*R*l2*sin(t2) - (l1*l2*cos(t1 - t2))/2 - (3*l1*l2*cos(t1 + t2))/2 - 3^(1/2)*a56*l1*sin(t1) + 3^(1/2)*a56*l2*sin(t2))^(1/2) - 3^(1/2)*r;(3*R^2 + 3*a56^2 + l2^2 + l3^2 - 3*R*l2*sin(t2) - 3*R*l3*sin(t3) - (l2*l3*cos(t2 - t3))/2 - (3*l2*l3*cos(t2 + t3))/2 - 3^(1/2)*a56*l2*sin(t2) + 3^(1/2)*a56*l3*sin(t3))^(1/2) - 3^(1/2)*r;(3*R^2 + 3*a56^2 + l1^2 + l3^2 - 3*R*l1*sin(t1) - 3*R*l3*sin(t3) - (l1*l3*cos(t1 - t3))/2 - (3*l1*l3*cos(t1 + t3))/2 + 3^(1/2)*a56*l1*sin(t1) - 3^(1/2)*a56*l3*sin(t3))^(1/2) - 3^(1/2)*r];
    thetas = thetas - jac\eq;
    iter = iter + 1;
end

theta1 = thetas(1); theta2 = thetas(2); theta3 = thetas(3);

P_c = [(l1*cos(theta1))/3 + (l2*cos(theta2))/3 + (l3*cos(theta3))/3;(R*cos(alpha_5))/3 + (a56*sin(alpha_5))/3 - (R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 - (R*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 - (l1*cos(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (l3*sin(theta3)*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3;(R*sin(alpha_5))/3 - (a56*cos(alpha_5))/3 - (R*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (l1*sin(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (l3*sin(theta3)*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3];

px = P_c(1);
py = P_c(2);
pz = P_c(3);

A_lin = [r*sin(alpha_13),0,0,r*cos(alpha_13),0,0;0,r*sin(alpha_13),0,0,r*cos(alpha_13),0;0,0,r*sin(alpha_13),0,0,r*cos(alpha_13);-r*(sin(alpha_13)/2+(3^(1/2)*cos(alpha_13))/2),0,0,-r*(cos(alpha_13)/2-(3^(1/2)*sin(alpha_13))/2),0,0;0,-r*(sin(alpha_13)/2+(3^(1/2)*cos(alpha_13))/2),0,0,-r*(cos(alpha_13)/2-(3^(1/2)*sin(alpha_13))/2),0;0,0,-r*(sin(alpha_13)/2+(3^(1/2)*cos(alpha_13))/2),0,0,-r*(cos(alpha_13)/2-(3^(1/2)*sin(alpha_13))/2)];
B_lin = [l1*cos(theta1) - px; R*cos(alpha_5) - py + a56*sin(alpha_5) - l1*cos(alpha_5)*sin(theta1); R*sin(alpha_5) - a56*cos(alpha_5) - pz - l1*sin(alpha_5)*sin(theta1); l2*cos(theta2) - px; l2*sin(theta2)*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2) - py; a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2) - R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2) - pz + l2*sin(theta2)*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2)];
X = linsolve(A_lin,B_lin);

X(1:3) = X(1:3)/norm(X(1:3));
X(4:6) = X(4:6)/norm(X(4:6));

a1_ = X(1);
a2_ = X(2);
a3_ = X(3);
o1_ = X(4);
o2_ = X(5);
o3_ = X(6);
n1_ =  o2_*a3_ - a2_*o3_;
n2_ = -o1_*a3_ + a1_*o3_;
n3_ =  o1_*a2_ - o2_*a1_;

% beta = asin(n2_); % CHECK THIS SIGN DO ATAN2S
gamma = atan2(-a2_,o2_);
% alpha = acos(n1_/cos(beta)); % CHECK THIS SIGN DO ATAN2S
alpha = atan2(-n3_,n1_);
% gamma = acos(o2_/cos(beta)); % CHECK THIS SIGN DO ATAN2S
beta = atan2(n2_,sqrt(n3_^2+n1_^2));

q_out = [pi/2-theta1 pi/2-theta2 pi/2-theta3 l1 l2 l3 px py pz alpha beta gamma].';
end