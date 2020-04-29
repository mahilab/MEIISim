%% GETTING W EQUATIONS FOR A SYMMETRIC PLATFORM
clear all
syms theta1 theta2 theta3 l1 l2 l3... 
     theta1_dot theta2_dot theta3_dot l1_dot l2_dot l3_dot...  
     px py pz vx vy vz wx wy wz...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     p         

% Vectors for generalized coordinates
l         = [        l1,         l2,         l3];
theta     = [    theta1,     theta2,     theta3];
l_dot     = [    l1_dot,     l2_dot,     l3_dot];
theta_dot = [theta1_dot, theta2_dot, theta3_dot];

% Velocity of platform
Vc = [vx; vy; vz];

% Angular velocity of platform
w  = [wx; wy; wz];

% Transformation matrix from ijk to IJK
T = [n1 o1 a1 px;
     n2 o2 a2 py;
     n3 o3 a3 pz;
      0  0  0  1];

Tinv = [T(1:3,1:3).',-T(1:3,1:3).'*T(1:3,4);T(4,:)];
% Eqn 27
r{1} = [   1;          0; 0]*p;
r{2} = [-1/2;  sqrt(3)/2; 0]*p;
r{3} = [-1/2; -sqrt(3)/2; 0]*p;
 
% Eqn 29
phi_l{1} = [   -cos(theta1);                      0; sin(theta1)];
phi_l{2} = [1/2*cos(theta2); -sqrt(3)/2*cos(theta2); sin(theta2)];
phi_l{3} = [1/2*cos(theta3);  sqrt(3)/2*cos(theta3); sin(theta3)];

% Eqn 30
phi_z{1} = [0;             1; 0];
phi_z{2} = [-sqrt(3)/2; -1/2; 0];
phi_z{3} = [ sqrt(3)/2; -1/2; 0];

for i = 1:3
    % Eqn 28
    Vbi_1{i} = l_dot(i)*phi_l{i} + cross(theta_dot(i)*phi_z{i},l(i)*phi_l{i});
    Temp{i} = T*[r{i};1];
    
    % Eqn 26
    Vbi_2{i} = Vc + cross(w,Temp{i}(1:3));
    
    % Equaliy between eqn 26 and eqn 28
    eq{i} = Vbi_1{i} == Vbi_2{i};
end

% Solving linear system of equations, choosing X relation for first joint,
% Y relation for second joint, and Z relation for third joint. Here, if I
% choose X,Y,Z relations for same joint, there is no solution.
[A,B] = equationsToMatrix([eq{1}(2);eq{3}(2)], [wx]);
w_solved = A\B;

wx_solved = simplify(w_solved(1))
wy_solved = simplify(w_solved(2))
% wz_solved = simplify(w_solved(3))

% wz_solved = solve(eq{1}(2),wz)
% wy_solved = solve(subs(eq{1}(1),wz,wz_solved),wy)
% wx_soln   = solve(subs(eq{1}(3),wy,wy_solved),wx)

%% Transforming Vc instead of Ro
clear all
syms theta1 theta2 theta3 l1 l2 l3... 
     theta1_dot theta2_dot theta3_dot l1_dot l2_dot l3_dot...  
     px py pz vx vy vz wx wy wz...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     p r        

gamma = [0, 2*pi/3, -2*pi/3];
 
% Vectors for generalized coordinates
l         = [        l1,         l2,         l3];
theta     = [    theta1,     theta2,     theta3];
l_dot     = [    l1_dot,     l2_dot,     l3_dot];
theta_dot = [theta1_dot, theta2_dot, theta3_dot];

% Velocity of platform
Vc = [vx; vy; vz];

% Angular velocity of platform
w  = [wx; wy; wz];

% Transformation matrix from ijk to IJK
T = [n1 o1 a1 px;
     n2 o2 a2 py;
     n3 o3 a3 pz;
      0  0  0  1];

Tinv = [T(1:3,1:3).',-T(1:3,1:3).'*T(1:3,4);T(4,:)];
% Eqn 27
r_{1} = [   1;          0; 0]*p;
r_{2} = [-1/2;  sqrt(3)/2; 0]*p;
r_{3} = [-1/2; -sqrt(3)/2; 0]*p;
 
% Eqn 29
phi_l{1} = [   -cos(theta1);                      0; sin(theta1)];
phi_l{2} = [1/2*cos(theta2); -sqrt(3)/2*cos(theta2); sin(theta2)];
phi_l{3} = [1/2*cos(theta3);  sqrt(3)/2*cos(theta3); sin(theta3)];

% Eqn 30
phi_z{1} = [0;             1; 0];
phi_z{2} = [-sqrt(3)/2; -1/2; 0];
phi_z{3} = [ sqrt(3)/2; -1/2; 0];

for i = 1:3
    r_temp{i} = Rz(gamma(i))*TRANSx(r);
    r_comp{i} = r_temp{i}(1:3,4);
    R_0_2{i} = Rz(gamma(i))*Ry(-pi/2-(pi/2-theta(i)));
    phi_l_comp{i} = R_0_2{i}(1:3,1);
    phi_z_comp{i} = R_0_2{i}(1:3,2);
end

for i = 1:3
    % Eqn 28
    Vbi_1{i} = Tinv(1:3,1:3)*(l_dot(i)*phi_l{i} + cross(theta_dot(i)*phi_z{i},l(i)*phi_l{i}));
%     Temp{i} = T*[r{i};1];
    
    % Eqn 26
    Vbi_2{i} = Tinv(1:3,1:3)*Vc + cross(w,r{i});
    
    % Equaliy between eqn 26 and eqn 28
    eq{i} = Vbi_1{i} == Vbi_2{i};
end

% Solving linear system of equations, choosing X relation for first joint,
% Y relation for second joint, and Z relation for third joint. Here, if I
% choose X,Y,Z relations for same joint, there is no solution.
wy_solved = solve(eq{1}(3),wy);
% wz_solved = simplify(w_solved(3))

% wz_solved = solve(eq{1}(2),wz)
% wy_solved = solve(subs(eq{1}(1),wz,wz_solved),wy)
% wx_soln   = solve(subs(eq{1}(3),wy,wy_solved),wx)

%%
syms px py pz
Pc = [px; py; pz];
n1 = (1 - l1*cos(theta1)-px)/p;
n2 = -py/p;
n3 = (l1*sin(theta1)-pz)/p;

o1 = n2;
o2 = (sqrt(3) - sqrt(3)*l2*cos(theta2)-3*py)/(sqrt(3)*p);
o3 = (2*l2*sin(theta2)+l1*sin(theta1)-3*pz)/(sqrt(3)*p);

a1 = n2*o3-o2*n3;
a2 = -n1*o3+o1*n3;
a3 = n1*o2-n2*o1;

unit_vecs = [n1 n2 n3 o1 o2 o3 a1 a2 a3];

variables = [p     px   py  pz    vx     vy    vz    l1   l2   l3 theta1 theta2 theta3    l1_dot    l2_dot     l3_dot    theta1_dot     theta2_dot theta3_dot];
values    = [0.5 0.01 0.01 0.1  0.02  -0.02  0.01  0.11 0.11 0.11   pi/4   pi/4   pi/4      0.01      0.02      -0.01          0.01          -0.02          0];

test1 = w_solved(2);
test2 = -1/p*(a1*(-l1_dot*cos(theta1) + l1*sin(theta1)*theta1_dot-vx) - a2*vy + a3*(l1_dot*sin(theta1)+l1*cos(theta1)*theta1_dot-vz));

vpa(subs(subs(test1),variables,values),4)
vpa(subs(subs(test2),variables,values),4)

% simplify(expand(subs(test1)))
% simplify(expand(subs(test2)))