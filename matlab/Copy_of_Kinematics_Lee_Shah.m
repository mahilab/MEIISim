clear all;

syms xc yc zc...
     y z...
     a56 alpha_5 alpha_13...
     r R...
     a_ b_ c_ gam

variables = [r R alpha_5 alpha_13 a56];
values = [0.05288174521 0.1044956 0.094516665 5*pi/180 -3.8340e-04];

gamma = [0, -2*pi/3, 2*pi/3];

T = Ry(a_)*Rz(b_)*Rx(c_);
T(1:3,4) = [xc,yc,zc];

for i = 1:3
    P_temp = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56);
    P{i} = P_temp(1:3,4);
    b_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    b{i} = b_temp{i}(1:3,4);
    B{i} = T*[b{i};1];
end

for i = 1:3
    eq1 = solve(y == P{i}(2),R);
    eq2 = solve(z == P{i}(3),R);
    z_y{i} = simplify(expand(solve(eq1 == eq2,z))); % Plane that each link is constrained to
end

for i = 1:3
    y_temp = T(2,1:3)*b{i}+yc;
    sys_eq{i} = T(3,1:3)*b{i} + zc == subs(z_y{i},[y],y_temp);
end

[A_lin,B_lin] = equationsToMatrix([sys_eq{2}, sys_eq{3}], [yc, zc]);

X = linsolve(A_lin,B_lin);

yc_solved = simplify(expand(X(1)))
zc_solved = simplify(expand(X(2)))
soln1 = subs(sys_eq{1},[yc,zc],[yc_solved,zc_solved]);

soln_c_ = simplify(expand(solve(simplify(expand(soln1)),c_,'Real',true)))

for i = 1:3
    q{i} = simplify(expand(sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2 + (B{i}(3) - P{i}(3))^2)));
    theta{i} = expand(simplify(atan(sqrt((B{i}(1) - P{i}(1))^2),sqrt((B{i}(2) - P{i}(2))^2 + (B{i}(3) - P{i}(3))^2))));
end

vpa(wrapToPi((subs(soln_c_(2),[a_,b_,variables],[0, 0,values]))),12)

%%
% a_ = 0;
% b_ = 0;
% c_ = 0;
% xc = 0.11985868974;
% yc = 0;
% zc = 0;
% 
% R = 0.1044956; % [m]
% r = 0.05288174521; % [m]
% a5 = 0.0268986; % [m]
% a6 = 0.0272820; % [m]
% a56 = -(a5-a6); % [m]
% alpha_5 = 0.094516665; % [rad]
% alpha_13 = 5*pi/180; % [rad]
syms a_ b_ c_ xc yc zc r R a56 alpha_5 alpha_13
variables = [a_ b_ c_ xc yc zc r R a56 alpha_5 alpha_13];
values = [0 0 0 0.11985873568 0 0 0.05288174521 0.1044956 3.8340e-04 0.094516665 5*pi/180 ];
% 0.11985873568mm
P1_test = vpa(subs(P{1},variables,values),32);
B1_test = vpa(subs(B{1}(1:3),variables,values),32);

P2_test = vpa(subs(P{2},variables,values),32);
B2_test = vpa(subs(B{2}(1:3),variables,values),32);

P3_test = vpa(subs(P{3},variables,values),32);
B3_test = vpa(subs(B{3}(1:3),variables,values),32);

% B1 = [0.11985873568;  0.05268051420;  0.00460894777];
% B2 = [0.11985873568; -0.02234879124; -0.04792713747];
% B3 = [0.11985873568; -0.03033172296;  0.04331818970];
% 
% P1 = [0; 0.10406538063; 0.00948018821];
% P2 = [0; -.04382260649; -.09486335739];
% P3 = [0; -0.06024277414; 0.08538316918];

l1 = sqrt((P1_test(1)-B1_test(1))^2 + (P1_test(2)-B1_test(2))^2 + (P1_test(3)-B1_test(3))^2)
% norm(P1_test-B1_test)
% l1 = sqrt((P1(1)-B1(1))^2 + (P1(2)-B1(2))^2 + (P1(3)-B1(3))^2)

% l1_test-.1305
% B1_test - B1
% B2_test - B2
% B3_test - B3
% P1_test - P1
% P2_test - P2
% P3_test - P3