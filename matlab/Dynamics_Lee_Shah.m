clear all;

syms l1 l2 l3 theta1 theta2 theta3...
     l1_dot l2_dot l3_dot theta1_dot theta2_dot theta3_dot...
     l1_d_dot l2_d_dot l3_d_dot theta1_d_dot theta2_d_dot theta3_d_dot...
     l1t(t) l2t(t) l3t(t) theta1t(t) theta2t(t) theta3t(t)...
     l1_dott(t) l2_dott(t) l3_dott(t) theta1_dott(t) theta2_dott(t) theta3_dott(t)...
     wx wy wz vx vy vz px py pz...
     r R alpha_5 a56 alpha_13...
     F1 F2 F3...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     g
 
T = [n1 o1 a1 px;
     n2 o2 a2 py;
     n3 o3 a3 pz;
      0  0  0  1];
 
Tinv = [T(1:3,1:3).',-T(1:3,1:3).'*T(1:3,4);T(4,:)];
 
gamma = [0, -2*pi/3, 2*pi/3];

qs       = [      l1,       l2,       l3,       theta1,       theta2,       theta3];
q_dots   = [  l1_dot,   l2_dot,   l3_dot,   theta1_dot,   theta2_dot,   theta3_dot];
q_d_dots = [l1_d_dot, l2_d_dot, l3_d_dot, theta1_d_dot, theta2_d_dot, theta3_d_dot];

qst     = [    l1t,     l2t,     l3t,     theta1t,     theta2t,     theta3t];
q_dotst = [l1_dott, l2_dott, l3_dott, theta1_dott, theta2_dott, theta3_dott];
 
ls =     [l1,         l2,     l3];
l_dots = [l1_dot, l2_dot, l3_dot];

thetas     = [    theta1,     theta2,     theta3];
theta_dots = [theta1_dot, theta2_dot, theta3_dot];

Q = [F1 F2 F3 0 0 0];

w = [wx wy wz].';

P_c = [(l1*cos(theta1))/3 + (l2*cos(theta2))/3 + (l3*cos(theta3))/3;(R*cos(alpha_5))/3 + (a56*sin(alpha_5))/3 - (R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 - (R*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 - (l1*cos(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (l3*sin(theta3)*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3;(R*sin(alpha_5))/3 - (a56*cos(alpha_5))/3 - (R*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (l1*sin(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (l3*sin(theta3)*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3];

px = P_c(1); py = P_c(2); pz = P_c(3);

V_c = [0 0 0].';

for i = 1:3
    V_c = V_c + diff(P_c,ls(i))*l_dots(i) + diff(P_c,thetas(i))*theta_dots(i);
end

vx = V_c(1); vy = V_c(2); vz = V_c(3);

for i = 1:3
    r_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    r_{i} = r_temp{i}(1:3,4);
    R_0_2{i} = Rx(alpha_5 + gamma(i))*Rz(-pi/2-(pi/2-thetas(i)));
    phi_l{i} = R_0_2{i}(1:3,1);
    phi_z{i} = R_0_2{i}(1:3,3);
end

for i = 1:3
   Vbi_1{i} = Tinv(1:3,1:3)*V_c + cross(w,r_{i});
   Vbi_2{i} = Tinv(1:3,1:3)*(l_dots(i)*phi_l{i} + cross(theta_dots(i)*phi_z{i},ls(i)*phi_l{i}));
   eq{i} = Vbi_1{i} - Vbi_2{i} == 0;
end

wx_solved = simplify(solve(eq{1}(2),wx));
[A,B] = equationsToMatrix([eq{1}(1);eq{2}(1)], [wy wz]);
wy_wz = A\B;
wy_solved = simplify(wy_wz(1));
wz_solved = simplify(wy_wz(2));

% substitue in unit vectors to get into a form we can use for dynamics
unit_vecs = [ -(abs(r)*(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1)))/(r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), (3^(1/2)*abs(r)*(3*a56*cos(alpha_5)*cos(alpha_13) - 3*py*sin(alpha_13) + 3*R*cos(alpha_5)*sin(alpha_13) - 3*R*cos(alpha_13)*sin(alpha_5) + 3*a56*sin(alpha_5)*sin(alpha_13) + 3*3^(1/2)*py*cos(alpha_13) - 3*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2)))/(3*r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*pz*sin(alpha_13) - 3*R*cos(alpha_5)*cos(alpha_13) + 3*a56*cos(alpha_5)*sin(alpha_13) - 3*a56*cos(alpha_13)*sin(alpha_5) - 3*R*sin(alpha_5)*sin(alpha_13) - 3*3^(1/2)*pz*cos(alpha_13) + 3*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2)))/(3*r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), (abs(r)*(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1)))/(r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*py*cos(alpha_13) + 3*3^(1/2)*py*sin(alpha_13) - 3*R*cos(alpha_5)*cos(alpha_13) + 3*a56*cos(alpha_5)*sin(alpha_13) - 3*a56*cos(alpha_13)*sin(alpha_5) - 3*R*sin(alpha_5)*sin(alpha_13) + 3*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2)))/(3*r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*pz*cos(alpha_13) + 3*3^(1/2)*pz*sin(alpha_13) + 3*a56*cos(alpha_5)*cos(alpha_13) + 3*R*cos(alpha_5)*sin(alpha_13) - 3*R*cos(alpha_13)*sin(alpha_5) + 3*a56*sin(alpha_5)*sin(alpha_13) + 3*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2)))/(3*r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*R^2 + 3*a56^2 - 3*R*l1*sin(theta1) - 3*R*l2*sin(theta2) - 3*a56*py*sin(alpha_5) - 3*R*py*cos(alpha_5) + 3*a56*pz*cos(alpha_5) - 3*R*pz*sin(alpha_5) + 3*l1*l2*sin(theta1)*sin(theta2) + 3*3^(1/2)*R*pz*cos(alpha_5) + 3*3^(1/2)*a56*py*cos(alpha_5) - 3*3^(1/2)*R*py*sin(alpha_5) + 3*3^(1/2)*a56*pz*sin(alpha_5) + 3^(1/2)*a56*l1*sin(theta1) - 3^(1/2)*a56*l2*sin(theta2) + 3*l2*py*cos(alpha_5)*sin(theta2) + 3*l2*pz*sin(alpha_5)*sin(theta2) - 2*3^(1/2)*l1*pz*cos(alpha_5)*sin(theta1) - 3^(1/2)*l2*pz*cos(alpha_5)*sin(theta2) + 2*3^(1/2)*l1*py*sin(alpha_5)*sin(theta1) + 3^(1/2)*l2*py*sin(alpha_5)*sin(theta2)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*a56*px*sin(alpha_5) + 3*R*px*cos(alpha_5) - 3*3^(1/2)*a56*px*cos(alpha_5) + 3*3^(1/2)*R*px*sin(alpha_5) - 2*3^(1/2)*l1*pz*cos(theta1) + 2*3^(1/2)*l2*pz*cos(theta2) - 3*R*l1*cos(alpha_5)*cos(theta1) - 3*a56*l1*sin(alpha_5)*cos(theta1) - 3*l2*px*cos(alpha_5)*sin(theta2) + 3^(1/2)*a56*l1*cos(alpha_5)*cos(theta1) + 2*3^(1/2)*a56*l2*cos(alpha_5)*cos(theta2) - 3^(1/2)*R*l1*sin(alpha_5)*cos(theta1) - 2*3^(1/2)*R*l2*sin(alpha_5)*cos(theta2) - 2*3^(1/2)*l1*px*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*px*sin(alpha_5)*sin(theta2) + 3*l1*l2*cos(alpha_5)*cos(theta1)*sin(theta2) + 3^(1/2)*l1*l2*sin(alpha_5)*cos(theta1)*sin(theta2) + 2*3^(1/2)*l1*l2*sin(alpha_5)*cos(theta2)*sin(theta1)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*R*px*sin(alpha_5) - 3*a56*px*cos(alpha_5) - 3*3^(1/2)*R*px*cos(alpha_5) - 3*3^(1/2)*a56*px*sin(alpha_5) + 2*3^(1/2)*l1*py*cos(theta1) - 2*3^(1/2)*l2*py*cos(theta2) + 3*a56*l1*cos(alpha_5)*cos(theta1) - 3*R*l1*sin(alpha_5)*cos(theta1) - 3*l2*px*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*l1*cos(alpha_5)*cos(theta1) + 2*3^(1/2)*R*l2*cos(alpha_5)*cos(theta2) + 3^(1/2)*a56*l1*sin(alpha_5)*cos(theta1) + 2*3^(1/2)*a56*l2*sin(alpha_5)*cos(theta2) + 2*3^(1/2)*l1*px*cos(alpha_5)*sin(theta1) + 3^(1/2)*l2*px*cos(alpha_5)*sin(theta2) + 3*l1*l2*sin(alpha_5)*cos(theta1)*sin(theta2) - 3^(1/2)*l1*l2*cos(alpha_5)*cos(theta1)*sin(theta2) - 2*3^(1/2)*l1*l2*cos(alpha_5)*cos(theta2)*sin(theta1)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2))];
unit_vecs = simplify(unit_vecs);
wx_solved = simplify(subs(wx_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
wy_solved = simplify(subs(wy_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
wz_solved = simplify(subs(wz_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));

%%
syms Ixx Iyy Izz Mp Ml l_offset;

sum_d_thetadot2  = 0;
sum_By  = 0;

for i = 1:3
    d{i} = (ls(i) - l_offset);
    sum_d_thetadot2 = sum_d_thetadot2 + d{i}^2*theta_dots(i)^2;
    B_temp{i} = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-thetas(i))*TRANSx(d{i});
    By{i} = B_temp{i}(2,4);
    sum_By = sum_By + By{i};
end

fk = [(abs(l1*cos(theta1) - l2*cos(theta2))^2 + abs((3*R*cos(alpha_5))/2 + (3*a56*sin(alpha_5))/2 - l1*cos(alpha_5)*sin(theta1) - (l2*cos(alpha_5)*sin(theta2))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l2*sin(alpha_5)*sin(theta2))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(theta1) + (l2*sin(alpha_5)*sin(theta2))/2 - (3^(1/2)*R*cos(alpha_5))/2 - (3^(1/2)*a56*sin(alpha_5))/2 + (3^(1/2)*l2*cos(alpha_5)*sin(theta2))/2)^2)^(1/2) - 3^(1/2)*r;(abs(l2*cos(theta2) - l3*cos(theta3))^2 + abs((l3*cos(alpha_5)*sin(theta3))/2 - (l2*cos(alpha_5)*sin(theta2))/2 + 3^(1/2)*a56*cos(alpha_5) - 3^(1/2)*R*sin(alpha_5) + (3^(1/2)*l2*sin(alpha_5)*sin(theta2))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(theta3))/2)^2 + abs((l2*sin(alpha_5)*sin(theta2))/2 - (l3*sin(alpha_5)*sin(theta3))/2 - 3^(1/2)*R*cos(alpha_5) - 3^(1/2)*a56*sin(alpha_5) + (3^(1/2)*l2*cos(alpha_5)*sin(theta2))/2 + (3^(1/2)*l3*cos(alpha_5)*sin(theta3))/2)^2)^(1/2) - 3^(1/2)*r; (abs(l1*cos(theta1) - l3*cos(theta3))^2 + abs(l1*cos(alpha_5)*sin(theta1) - (3*a56*sin(alpha_5))/2 - (3*R*cos(alpha_5))/2 + (l3*cos(alpha_5)*sin(theta3))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(theta3))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(theta1) + (l3*sin(alpha_5)*sin(theta3))/2 + (3^(1/2)*R*cos(alpha_5))/2 + (3^(1/2)*a56*sin(alpha_5))/2 - (3^(1/2)*l3*cos(alpha_5)*sin(theta3))/2)^2)^(1/2) - 3^(1/2)*r];

T = 1/2*Mp*(vx^2 + vy^2 + vz^2) + 1/2*(Ixx*wx_solved^2 + Iyy*wy_solved^2 + Izz*wz_solved^2) + 1/2*Ml*sum_d_thetadot2;

P = Mp*g*py + simplify(Ml*g*sum_By);

L = T-P;
L = subs(L,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs);

%% COMPUTING EOMS
syms gamm_1 gamm_2 gamm_3

for i = 1:6
    dldqdot{i} = diff(L,q_dots(i));
    dldqdot_t{i} = subs(dldqdot{i},[qs q_dots], [qst q_dotst]);
    d_dldqdot_t{i} = diff(dldqdot_t{i},t);
    d_dldqdot{i} = subs(d_dldqdot_t{i},[qst q_dotst diff([qst q_dotst],t)],[qs q_dots q_dots q_d_dots]);
    
    EOMS(i) = d_dldqdot{i} - diff(L,qs(i));% + gamm_1*diff(fk(1),qs(i)) + gamm_2*diff(fk(2),qs(i)) + gamm_3*diff(fk(3),qs(i));
end

[M V G] = separate_mvg_no_simp(EOMS.',q_d_dots.',g);
%%
Vsquare = sym(zeros(6,6));

for i = 1:6
    test = V(i);
    for j = 1:6
        [coe,terms] = coeffs(test,q_dots(j));
        for k = 1:length(terms)
            if subs(terms(k),q_dots(j),0) ~= 0
                Vsquare(i,j) = Vsquare(i,j) + coe(k)*terms(k);
            else
                test = coe(k);
            end
        end
    end
end

%%

psi = [fk;qs(1:3).'];
psi_dq = simplify(jacobian(psi,qs));
rho = simplify(inv(psi_dq)*[zeros(3,3);eye(3,3)]);
rho_t = subs(rho,qs,qst);
rho_dt = subs(diff(rho_t,t),[qst diff(qst,t)],[qs, q_dots]);

Mpar = rho.'*M*rho;
Vpar = rho.'*Vsquare*rho + rho.'*M*rho_dt;
Gpar = rho.'*G;


%%
% AT = simplify(jacobian(fk,[theta1 theta2 theta3]));
% AinvT = simplify(inv(AT));
% BT = simplify(jacobian(fk,[l1 l2 l3]).');
% gammas = AinvT*EOMS(4:6).';
% EOMS_final = [EOMS(1:3).' - BT*gammas; EOMS(4:6).' - AT*gammas];
%% Substitute in actual values

% Ixx = 0.03990366;
% Iyy = 0.09672265;
% Izz = 0.13629088;
Ixx = 0.00091827;
Iyy = 0.00080377;
Izz = 0.00138661;
Mp = 0.36486131;
Ml = 0.14820004;
l_offset = 0.08803884;
r        = 0.05288174521;
R        = 0.1044956;
alpha_5  = 0.094516665; 
alpha_13 = 5*pi/180;
a56      = 3.8340e-04;

Mpar = subs(Mpar);
Vpar = subs(Vpar);
Gpar = subs(Gpar);
% g        = 9.81;

% EOMS_vals = vpa(subs(EOMS_final),5);

%% Convert to qdd form
% [M V G] = separate_mvg_no_simp(EOMS_final,q_d_dots.',g);

g = 9.80665;

Mpar = subs(Mpar);
Vpar = subs(Vpar);
Gpar = subs(Gpar);

%% Save to mat file instead
% M = vpa(M,3);
% V = vpa(V,3);
% G = vpa(G,3);
% 
% save('MVG','M','V','G')

%% Testing
% mass_props = [Ixx Iyy Izz Mp Ml l_offset g];
% mass_values = [0.03990366 0.09672265 0.13629088 0.36486131 0.14820004 0.08803884 9.81];
% q_q_dot_vals = [0.1 0.1 0.11 1.0390 0.9300 1.0663 0.001 0.0001 -0.0003 0.001 0.001 -0.003];
% variables = [r R alpha_5 alpha_13 a56];
% values = [0.05288174521 0.1044956 0.094516665 5*pi/180 3.8340e-04];
% M_eval = subs(subs(M),[qs q_dots variables mass_props], [q_q_dot_vals values mass_values]);
% vpa(M_eval,2)
% V_eval = subs(subs(V),[qs q_dots variables mass_props], [q_q_dot_vals values mass_values]);
% vpa(V_eval,2)
% G_eval = subs(subs(G),[qs q_dots variables mass_props], [q_q_dot_vals values mass_values]);
% vpa(G_eval,2)

%% Write to file

% for i = 1:6
%     for j = 1:6
%         M_write = char(M(i,j),3);
%         fileID = fopen("DynamicEqs/M" + int2str(i) + int2str(j) + ".txt","w");
%         fprintf(fileID,M_write);
%     end
%     V_write = char(V(i));
%     fileID = fopen("DynamicEqs/V" + int2str(i) + ".txt","w");
%     fprintf(fileID,V_write);
%     G_write = char(G(i));
%     fileID = fopen("DynamicEqs/G" + int2str(i) + ".txt","w");
%     fprintf(fileID,G_write);
% end

%% Write parallel structure to file

for i = 1:3
    for j = 1:3
        M_write = ccode(vpa(Mpar(i,j),3));
        fileID = fopen("DynamicEqs/Mpar" + int2str(i) + int2str(j) + ".txt","w");
        fprintf(fileID,M_write);
        V_write = ccode(vpa(Vpar(i,j),3));
        fileID = fopen("DynamicEqs/Vpar" + int2str(i) + int2str(j) + ".txt","w");
        fprintf(fileID,V_write);
    end
    G_write = ccode(vpa(Gpar(i),3));
    fileID = fopen("DynamicEqs/Gpar" + int2str(i) + ".txt","w");
    fprintf(fileID,G_write);
end
%%
% mass_props = [Ixx Iyy Izz Mp Ml l_offset g];
% mass_values = [0.03990366 0.09672265 0.13629088 0.36486131 0.14820004 0.08803884 9.81];
% q_q_dot_vals = [0.1 0.1 0.11 1.0390 0.9300 1.0663 0.001 0.0001 -0.0003 0.001 0.001 -0.003];
% q_ddot_vals = [0.02 0.01 -0.01 0.03 -0.002 0.04];
% variables = [r R alpha_5 alpha_13 a56];
% values = [0.05288174521 0.1044956 0.094516665 5*pi/180 3.8340e-04];
% test = subs(subs(Vpar(1,1)),[qs q_dots q_d_dots variables mass_props], [q_q_dot_vals q_ddot_vals values mass_values]);
% 
% vpa(test,4)