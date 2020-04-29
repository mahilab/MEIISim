clear all;

syms l1(t) l2(t) l3(t) theta1(t) theta2(t) theta3(t) l1_dot(t) l2_dot(t) l3_dot(t) theta1_dot(t) theta2_dot(t) theta3_dot(t)...
     wx wy wz vx vy vz...
     r R alpha_5 a56 alpha_13...
     F1 F2 F3...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     g
 
T = [n1 o1 a1;
     n2 o2 a2;
     n3 o3 a3];
 
gamma = [0, -2*pi/3, 2*pi/3];

ls =     [        l1,         l2,         l3];
l_dots = [diff(l1,t), diff(l2,t), diff(l3,t)];

ls = ls(t);
l_dots = l_dots(t);

thetas     = [         theta1,        theta2,         theta3];
theta_dots = [diff(theta1,t), diff(theta2,t), diff(theta3,t)];

thetas = thetas(t);
theta_dots = theta_dots(t);

qs     = [       l1,        l2,        l3,        theta1,        theta2,        theta3];
q_dots = [l_dots(1), l_dots(2), l_dots(3), theta_dots(1), theta_dots(2), theta_dots(3)];

qs = qs(t);

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
    R_0_2{i} = Rx(alpha_5 + gamma(i))*Rz(-thetas(i));
    phi_l{i} = R_0_2{i}(1:3,1);
    phi_z{i} = R_0_2{i}(1:3,3);
end

% V_c = [vx vy vz].';

for i = 1:3
   Vbi_1{i} = V_c + cross(w,T*r_{i});
   Vbi_2{i} = l_dots(i)*phi_l{i} + cross(theta_dots(i)*phi_z{i},ls(i)*phi_l{i});
   eq{i} = Vbi_1{i} - Vbi_2{i} == 0;
end

% wz_solved = solve(eq{2}(1),wz);
% wx_solved = solve(subs(eq{2}(2),wz,wz_solved),wx);
% wy_soln = solve(subs(eq{1}(3),wx,wx_solved),wy)

[A,B] = equationsToMatrix([eq{2}(1);eq{1}(2);eq{1}(3)], [wx wy wz]);
w = A\B;
% unit_vecs = [ -(abs(r)*(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1)))/(r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), (3^(1/2)*abs(r)*(3*a56*cos(alpha_5)*cos(alpha_13) - 3*py*sin(alpha_13) + 3*R*cos(alpha_5)*sin(alpha_13) - 3*R*cos(alpha_13)*sin(alpha_5) + 3*a56*sin(alpha_5)*sin(alpha_13) + 3*3^(1/2)*py*cos(alpha_13) - 3*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2)))/(3*r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*pz*sin(alpha_13) - 3*R*cos(alpha_5)*cos(alpha_13) + 3*a56*cos(alpha_5)*sin(alpha_13) - 3*a56*cos(alpha_13)*sin(alpha_5) - 3*R*sin(alpha_5)*sin(alpha_13) - 3*3^(1/2)*pz*cos(alpha_13) + 3*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2)))/(3*r*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)), (abs(r)*(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1)))/(r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*py*cos(alpha_13) + 3*3^(1/2)*py*sin(alpha_13) - 3*R*cos(alpha_5)*cos(alpha_13) + 3*a56*cos(alpha_5)*sin(alpha_13) - 3*a56*cos(alpha_13)*sin(alpha_5) - 3*R*sin(alpha_5)*sin(alpha_13) + 3*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2)))/(3*r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), -(3^(1/2)*abs(r)*(3*pz*cos(alpha_13) + 3*3^(1/2)*pz*sin(alpha_13) + 3*a56*cos(alpha_5)*cos(alpha_13) + 3*R*cos(alpha_5)*sin(alpha_13) - 3*R*cos(alpha_13)*sin(alpha_5) + 3*a56*sin(alpha_5)*sin(alpha_13) + 3*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2)))/(3*r*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*R^2 + 3*a56^2 - 3*R*l1*sin(theta1) - 3*R*l2*sin(theta2) - 3*a56*py*sin(alpha_5) - 3*R*py*cos(alpha_5) + 3*a56*pz*cos(alpha_5) - 3*R*pz*sin(alpha_5) + 3*l1*l2*sin(theta1)*sin(theta2) + 3*3^(1/2)*R*pz*cos(alpha_5) + 3*3^(1/2)*a56*py*cos(alpha_5) - 3*3^(1/2)*R*py*sin(alpha_5) + 3*3^(1/2)*a56*pz*sin(alpha_5) + 3^(1/2)*a56*l1*sin(theta1) - 3^(1/2)*a56*l2*sin(theta2) + 3*l2*py*cos(alpha_5)*sin(theta2) + 3*l2*pz*sin(alpha_5)*sin(theta2) - 2*3^(1/2)*l1*pz*cos(alpha_5)*sin(theta1) - 3^(1/2)*l2*pz*cos(alpha_5)*sin(theta2) + 2*3^(1/2)*l1*py*sin(alpha_5)*sin(theta1) + 3^(1/2)*l2*py*sin(alpha_5)*sin(theta2)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*a56*px*sin(alpha_5) + 3*R*px*cos(alpha_5) - 3*3^(1/2)*a56*px*cos(alpha_5) + 3*3^(1/2)*R*px*sin(alpha_5) - 2*3^(1/2)*l1*pz*cos(theta1) + 2*3^(1/2)*l2*pz*cos(theta2) - 3*R*l1*cos(alpha_5)*cos(theta1) - 3*a56*l1*sin(alpha_5)*cos(theta1) - 3*l2*px*cos(alpha_5)*sin(theta2) + 3^(1/2)*a56*l1*cos(alpha_5)*cos(theta1) + 2*3^(1/2)*a56*l2*cos(alpha_5)*cos(theta2) - 3^(1/2)*R*l1*sin(alpha_5)*cos(theta1) - 2*3^(1/2)*R*l2*sin(alpha_5)*cos(theta2) - 2*3^(1/2)*l1*px*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*px*sin(alpha_5)*sin(theta2) + 3*l1*l2*cos(alpha_5)*cos(theta1)*sin(theta2) + 3^(1/2)*l1*l2*sin(alpha_5)*cos(theta1)*sin(theta2) + 2*3^(1/2)*l1*l2*sin(alpha_5)*cos(theta2)*sin(theta1)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2)), (abs(r)^2*(3*R*px*sin(alpha_5) - 3*a56*px*cos(alpha_5) - 3*3^(1/2)*R*px*cos(alpha_5) - 3*3^(1/2)*a56*px*sin(alpha_5) + 2*3^(1/2)*l1*py*cos(theta1) - 2*3^(1/2)*l2*py*cos(theta2) + 3*a56*l1*cos(alpha_5)*cos(theta1) - 3*R*l1*sin(alpha_5)*cos(theta1) - 3*l2*px*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*l1*cos(alpha_5)*cos(theta1) + 2*3^(1/2)*R*l2*cos(alpha_5)*cos(theta2) + 3^(1/2)*a56*l1*sin(alpha_5)*cos(theta1) + 2*3^(1/2)*a56*l2*sin(alpha_5)*cos(theta2) + 2*3^(1/2)*l1*px*cos(alpha_5)*sin(theta1) + 3^(1/2)*l2*px*cos(alpha_5)*sin(theta2) + 3*l1*l2*sin(alpha_5)*cos(theta1)*sin(theta2) - 3^(1/2)*l1*l2*cos(alpha_5)*cos(theta1)*sin(theta2) - 2*3^(1/2)*l1*l2*cos(alpha_5)*cos(theta2)*sin(theta1)))/(r^2*(abs(3^(1/2)*px*sin(alpha_13) - 3*px*cos(alpha_13) + l1*cos(alpha_13)*cos(theta1) + 2*l2*cos(alpha_13)*cos(theta2) - 3^(1/2)*l1*sin(alpha_13)*cos(theta1))^2 + abs(3*pz*cos(alpha_13) - 3^(1/2)*pz*sin(alpha_13) + l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - l2*cos(alpha_13)*sin(alpha_5)*sin(theta2) + 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) - 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) + 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) - 3^(1/2)*l1*sin(alpha_5)*sin(alpha_13)*sin(theta1))^2 + abs(3*py*cos(alpha_13) - 3^(1/2)*py*sin(alpha_13) + l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) - l2*cos(alpha_5)*cos(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) + 3^(1/2)*l2*cos(alpha_13)*sin(alpha_5)*sin(theta2))^2)^(1/2)*(abs(l1*sin(alpha_13)*cos(theta1) - 3*px*sin(alpha_13) + 2*l2*sin(alpha_13)*cos(theta2) - 3^(1/2)*px*cos(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*cos(theta1))^2 + abs(3*py*sin(alpha_13) + 3^(1/2)*py*cos(alpha_13) + l1*cos(alpha_5)*sin(alpha_13)*sin(theta1) - l2*cos(alpha_5)*sin(alpha_13)*sin(theta2) - 3^(1/2)*R*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*a56*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*a56*cos(alpha_13)*sin(alpha_5) - 3^(1/2)*R*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_5)*cos(alpha_13)*sin(theta1) + 3^(1/2)*l2*sin(alpha_5)*sin(alpha_13)*sin(theta2))^2 + abs(3*pz*sin(alpha_13) + 3^(1/2)*pz*cos(alpha_13) + l1*sin(alpha_5)*sin(alpha_13)*sin(theta1) - l2*sin(alpha_5)*sin(alpha_13)*sin(theta2) + 3^(1/2)*a56*cos(alpha_5)*cos(alpha_13) + 3^(1/2)*R*cos(alpha_5)*sin(alpha_13) - 3^(1/2)*R*cos(alpha_13)*sin(alpha_5) + 3^(1/2)*a56*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*l1*cos(alpha_13)*sin(alpha_5)*sin(theta1) - 3^(1/2)*l2*cos(alpha_5)*sin(alpha_13)*sin(theta2))^2)^(1/2))];
% w = subs(w,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs);
wx = simplify(w(1)); wy = simplify(w(2)); wz = simplify(w(3));

%%
syms Ixx Iyy Izz Mp Ml;

% Mp = 0.36486131; % kg

% Ixx = 0.03990366; % kg*m^2
% Iyy = 0.09672265; % kg*m^2
% Izz = 0.13629088; % kg*m^2

% Ml = 0.14820004; % kg

sum_d_thetadot2  = 0;
sum_By  = 0;

for i = 1:3
    d{i} = (ls(i) - 0.08803884);
    sum_d_thetadot2 = sum_d_thetadot2 + d{i}^2*theta_dots(i)^2;
    B_temp{i} = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-thetas(i))*TRANSx(d{i});
    By{i} = B_temp{i}(2,4);
    sum_By = sum_By + By{i};
end

T = 1/2*Mp*(vx^2 + vy^2 + vz^2) + 1/2*(Ixx*wx^2 + Iyy*wy^2 + Izz*wz^2) + 1/2*Ml*sum_d_thetadot2;

P = Mp*g*py + Ml*g*sum_By;

L = T-P;

%% COMPUTING EOMS

for i = 1:6
    
end
