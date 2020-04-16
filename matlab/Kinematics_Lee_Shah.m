%% Solution for symmetric 3rps mechanism based on 
% Kinematic Analysis of a Three-Degrees-of-Freedom In-Parallel Actuated
% Manipulator from the following source
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=796

clear all; clc;
syms r R...
     xc yc zc...
     n1 n2 n3...
     o1 o2 o3...
     a1 a2 a3...
     x y z...
     a_ b_ c_
 
gamma = [0, 2*pi/3, 4*pi/3];
 
% Eqn 3
T = Rz(a_)*Ry(b_)*Rz(c_);
T(1:3,4) = [xc,yc,zc];

for i = 1:3
    %Eqn 1
    P{i} = [cos(gamma(i))*R, sin(gamma(i))*R, 0].';
    %Eqn 2
    b{i} = [cos(gamma(i))*r, sin(gamma(i))*r, 0].';
    %Eqn 5
    B{i} = T*[b{i};1];
end

% Preceding eqn 9, this does not solve bc it is equal to 0
y_x{1} = 0;

for i = 2:3
    
    eq1 = solve(x == P{i}(1),R);
    eq2 = solve(y == P{i}(2),R);
    y_x{i} = simplify(expand(solve(eq1 == eq2,y)));
end

for i = 1:3
    x_temp = T(1,1:3)*b{i}+xc;
    sys_eq{i} = T(2,1:3)*b{i} + yc == subs(y_x{i},[x],x_temp);
end

[A_lin,B_lin] = equationsToMatrix([sys_eq{2}, sys_eq{3}], [xc, yc]);

X = linsolve(A_lin,B_lin);

xc_solved = simplify(expand(X(1)))
yc_solved = simplify(expand(X(2)))
soln1 = subs(sys_eq{1},[xc,yc],[xc_solved,yc_solved]);

soln_c_ = simplify(solve(simplify(expand(soln1)),c_,'Real',true))

vpa(subs(soln_c_(1),[a_,b_],[0.05,0.08]),6)

for i = 1:3
    q{i} = simplify(expand(sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2 + (B{i}(3) - P{i}(3))^2)));
    theta{i} = expand(simplify(atan(B{i}(3) - P{i}(3),sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2))));
end
%%
clear all;

syms xc yc zc...
     x y ...
     a56 alpha_5 alpha_13...
     r R...
     a_ b_ c_

variables = [r R alpha_5 alpha_13 a56];
values = [0.05288174521 0.1044956 0.094516665 5*pi/180 3.8340e-04];

gamma = [0, 2*pi/3, 4*pi/3];

T = Rx(a_)*Ry(b_)*Rz(c_);
T(1:3,4) = [xc,yc,zc];

for i = 1:3
    P_temp = Rz(alpha_5 + gamma(i))*TRANSx(R)*TRANSy(a56);
    P{i} = P_temp(1:3,4);
    b_temp = Rz(alpha_13 + gamma(i))*TRANSx(r);
    b{i} = b_temp(1:3,4);
    B{i} = T*[b{i};1];
end

for i = 1:3
    eq1 = solve(x == P{i}(1),R);
    eq2 = solve(y == P{i}(2),R);
    y_x{i} = simplify(expand(solve(eq1 == eq2,y))); % Plane that each link is constrained to
end

for i = 1:3
    x_temp = T(1,1:3)*b{i}+xc;
    sys_eq{i} = T(2,1:3)*b{i} + yc == subs(y_x{i},[x],x_temp);
end

[A_lin,B_lin] = equationsToMatrix([sys_eq{2}, sys_eq{3}], [xc, yc]);

X = linsolve(A_lin,B_lin);

xc_solved = simplify(expand(X(1)))
yc_solved = simplify(expand(X(2)))
soln1 = subs(sys_eq{1},[xc,yc],[xc_solved,yc_solved]);

soln_c_ = simplify(expand(solve(simplify(expand(soln1)),c_,'Real',true)))

for i = 1:3
    q{i} = simplify(expand(sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2 + (B{i}(3) - P{i}(3))^2)));
    theta{i} = expand(simplify(atan(B{i}(3) - P{i}(3),sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2))));
end

vpa(wrapToPi((subs(soln_c_(2),[a_,b_,variables],[-0.455173816243597, -0.469953566203639,values]))-pi),12)
%%
theta_0 = pi/4;
l_0 = 0.1305;
max_iter = 10;
tol = 1e-50;
par_sel = [4,5,6];
ser_sel = [10,11,7];
qp_0 = [theta_0; theta_0; theta_0; l_0; l_0; l_0; l_0*sin(theta_0); 0; 0; 0; 0; 0 ];

% q_ser = [-0.455173816243597, -0.469953566203639, 0.078685071468735 ]';
alpha = -25:25;
beta = -25:25;
x = 0.08:0.001:0.12;

diffs = zeros(12,size(alpha,2)*size(beta,2)*size(x,2));
i = 1;
for a = alpha
    for b = beta
        for xc = x
            q_ser = [deg2rad(a), deg2rad(b), xc];
            qp_ik = rps_solve(qp_0,ser_sel,q_ser,max_iter,tol);
            qp_ik_2 = MEII_IK_3(q_ser);
            diffs(1:12,i) = abs(qp_ik-qp_ik_2');
            i = i + 1;
        end
    end
end

mean_diffs = mean(diffs,2)
max_diffs = max(diffs,[],2)

%%
clear all;
theta_0 = pi/4;
l_0 = 0.11;
max_iter = 10;
tol = 1e-12;
par_sel = [4,5,6];
ser_sel = [10,11,7];
qp_0 = [theta_0; theta_0; theta_0; l_0; l_0; l_0; l_0*sin(theta_0); 0; 0; 0; 0; 0 ];

% q_ser = [-0.455173816243597, -0.469953566203639, 0.078685071468735 ]';
alpha = 0;
beta = 0;
x = 0.11985873568;

diffs = zeros(12,size(alpha,2)*size(beta,2)*size(x,2));
i = 1;
for a = alpha
    for b = beta
        for xc = x
            q_ser = [deg2rad(a), deg2rad(b), xc];
            qp_ik = rps_solve(qp_0,ser_sel,q_ser,max_iter,tol);
            qp_ik_2 = MEII_IK_3(q_ser);
            diffs(1:12,i) = abs(qp_ik-qp_ik_2');
            i = i + 1;
        end
    end
end

mean_diffs = mean(diffs,2)
% max_diffs = max(diffs,[],2)

l_corr = 0.13050000000; % m
theta_corr = 2.73495967-pi/2; % rad
abs(qp_ik(1)-theta_corr)
abs(qp_ik(4)-l_corr)

vpa(abs(qp_ik_2(1)-theta_corr),12)
vpa(abs(qp_ik_2(4)-l_corr),12)
