syms a b p Xc Yc Zc

c = -a;

euler = Rz(a)*Ry(b)*Rz(c);

n1 = euler(1,1);
n2 = euler(2,1);
n3 = euler(3,1);

o1 = euler(1,2);
o2 = euler(2,2);
o3 = euler(3,2);

a1 = euler(1,3);
a2 = euler(2,3);
a3 = euler(3,3);

L1_2 = (n1*p + Xc - 1)^2 + (n2*p+Yc)^2 + (n3*p+Zc)^2;

comp = 1 + p^2 + Xc^2 + Yc^2 + Zc^2 - 2*Xc + 2*p*(cos(a)^2*cos(b)+sin(a)^2)*(Xc-1) + p*(cos(b) - 1)*sin(2*a)*Yc - 2*p*sin(b)*cos(a)*Zc;

%% 
clear all; clc;
syms r R...
     xc yc zc...
     n1 n2 n3...
     o1 o2 o3...
     a1 a2 a3...
     x y z...
     a_ b_ c_
 
gamma = [0, 2*pi/3, 4*pi/3];

Xc = xc/R; Yc = yc/R; Zc = zc/R;
p = r/R;

% T = [n1 o1 a1 xc;
%      n2 o2 a2 yc;
%      n3 o3 a3 zc;
%      0  0  0  1];
 
T = Rz(a_)*Ry(b_)*Rz(c_);
T(1:3,4) = [xc,yc,zc];

for i = 1:3
    P{i} = [cos(gamma(i))*R, sin(gamma(i))*R, 0].';
    b{i} = [cos(gamma(i))*r, sin(gamma(i))*r, 0].';
    B{i} = T*[b{i};1];
end

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

[A,B] = equationsToMatrix([sys_eq{2}, sys_eq{3}], [xc, yc]);

X = linsolve(A,B);

xc_solved = simplify(expand(X(1)));
yc_solved = simplify(expand(X(2)));
soln1 = subs(sys_eq{1},[xc,yc],[xc_solved,yc_solved]);

soln_c_ = simplify(solve(simplify(expand(soln1)),c_,'Real',true))

vpa(subs(soln_c_(1),[a_,b_],[0.05,0.08]),6)

simplify(expand(subs(xc_solved,c_,soln_c_)))
simplify(expand(subs(yc_solved,c_,soln_c_)))

%%
clear all;

syms xc yc zc...
     x y ...
     a56 alpha_5 alpha_13...
     r R...
     a_ b_ c_
%      n1 n2 n3...
%      o1 o2 o3...
%      a1 a2 a3...

     
% alpha_13 = 5*pi/180;
% r        = 0.05288174521;
% R        = 0.1044956;   % meters
% alpha_5  = 0.094516665; % radians
% a4       = 0.159385;    % meters
% a5       = 0.0268986;   % meters
% a6       = 0.027282;    % meters
% a56      = 3.8340e-04;  % meters
variables = [r R alpha_5 alpha_13 a56];
values = [0.05288174521 0.1044956 0.094516665 5*pi/180 3.8340e-04];

gamma = [0, 2*pi/3, 4*pi/3];

% T = [n1 o1 a1 xc;
%      n2 o2 a2 yc;
%      n3 o3 a3 zc;
%      0  0  0  1];

T = Rx(a_)*Ry(b_)*Rz(c_);
T(1:3,4) = [xc,yc,zc];

% n1 = euler(1,1);
% n2 = euler(2,1);
% n3 = euler(3,1);
% 
% o1 = euler(1,2);
% o2 = euler(2,2);
% o3 = euler(3,2);
% 
% a1 = euler(1,3);
% a2 = euler(2,3);
% a3 = euler(3,3);

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

xc_solved = subs(simplify(expand(X(1))),variables,values)
yc_solved = subs(simplify(expand(X(2))),variables,values)
soln1 = vpa(simplify(expand(subs(sys_eq{1},[xc,yc,variables],[xc_solved,yc_solved,values]))),2);

soln_c_ = simplify(expand(solve(simplify(expand(soln1)),c_,'Real',true)))

vpa((subs(soln_c_(2),[a_,b_],[-0.455173816243597, -0.469953566203639])),12)
% L1_2 = ((P{1}(1) - B{1}(1))^2 + (P{1}(2) - B{1}(2))^2 + (P{1}(3) - B{1}(3))^2)/R^2

for i = 1:3
    q{i} = simplify(expand(subs(sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2 + (B{i}(3) - P{i}(3))^2),variables,values)));
    theta{i} = expand(simplify(subs(atan(B{i}(3) - P{i}(3),sqrt((B{i}(1) - P{i}(1))^2 + (B{i}(2) - P{i}(2))^2)),variables,values)));
end
%%
clear all;
sign1 = 1; sign2 = 1;
a_ = -0.455173816243597;
b_ = -0.469953566203639;

c_ = -2.0*atan((465140524608.0*cos(a_) + 465140524608.0*cos(b_) + 3372422040.0*sin(a_)*sin(b_) + (216367080863021331315264.0*cos(a_)^2*cos(b_)^2 + 432734161726042662630528.0*cos(a_)*cos(b_) + 216321587940696825546623.0)^(1/2))/(3372422040.0*cos(a_) + 3372422040.0*cos(b_) - 465140524608.0*sin(a_)*sin(b_) - 6744844129.0))
xc_solved = (3536240838930121*sin(212832608637151/2251799813685248))/4611686018427387904 - (1428949248134745*sin(pi/36)*cos(a_)*sin(c_))/18014398509481984 - (476316416044915*sin(pi/36)*cos(b_)*sin(c_))/18014398509481984 + (1428949248134745*cos(pi/36)*cos(a_)*cos(c_))/18014398509481984 + (476316416044915*cos(pi/36)*cos(b_)*cos(c_))/18014398509481984 - (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)^2*cos(a_)*cos(c_))/4503599627370496 + (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)^2*cos(a_)*sin(c_))/4503599627370496 - (1428949248134745*cos(pi/36)*sin(a_)*sin(b_)*sin(c_))/18014398509481984 - (1428949248134745*sin(pi/36)*cos(c_)*sin(a_)*sin(b_))/18014398509481984 + (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)^2*sin(a_)*sin(b_)*sin(c_))/4503599627370496 + (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)^2*cos(c_)*sin(a_)*sin(b_))/4503599627370496 - (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*cos(b_)*sin(c_))/4503599627370496 - (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*cos(b_)*cos(c_))/4503599627370496
yc_solved = (476316416044915*cos(pi/36)*cos(a_)*sin(c_))/18014398509481984 - (3536240838930121*cos(212832608637151/2251799813685248))/4611686018427387904 + (476316416044915*sin(pi/36)*cos(a_)*cos(c_))/18014398509481984 - (476316416044915*cos(pi/36)*cos(b_)*sin(c_))/18014398509481984 - (476316416044915*sin(pi/36)*cos(b_)*cos(c_))/18014398509481984 - (476316416044915*sin(pi/36)*sin(a_)*sin(b_)*sin(c_))/18014398509481984 + (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)^2*cos(b_)*sin(c_))/4503599627370496 + (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)^2*cos(b_)*cos(c_))/4503599627370496 + (476316416044915*cos(pi/36)*cos(c_)*sin(a_)*sin(b_))/18014398509481984 - (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*cos(a_)*cos(c_))/4503599627370496 + (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*cos(a_)*sin(c_))/4503599627370496 + (476316416044915*cos(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*sin(a_)*sin(b_)*sin(c_))/4503599627370496 + (476316416044915*sin(pi/36)*cos(212832608637151/2251799813685248)*sin(212832608637151/2251799813685248)*cos(c_)*sin(a_)*sin(b_))/4503599627370496

vpa(c_ - pi,12)
vpa(-1*xc_solved,12)
vpa(-1*yc_solved,12)

% q1 = (R^2 + a56^2 + r^2 + xc^2 + yc^2 + zc^2 + 2*a56*xc*sin(alpha_5) - 2*R*xc*cos(alpha_5) - 2*a56*yc*cos(alpha_5) - 2*R*yc*sin(alpha_5) + 2*r*xc*cos(alpha_13)*cos(b_)*cos(c_) + 2*r*yc*cos(a_)*cos(alpha_13)*sin(c_) + 2*r*yc*cos(a_)*cos(c_)*sin(alpha_13) - 2*r*xc*cos(b_)*sin(alpha_13)*sin(c_) + 2*r*zc*cos(alpha_13)*sin(a_)*sin(c_) + 2*r*zc*cos(c_)*sin(a_)*sin(alpha_13) - 2*R*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*cos(c_) - 2*a56*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*sin(c_) - 2*a56*r*cos(a_)*cos(alpha_5)*cos(c_)*sin(alpha_13) - 2*R*r*cos(a_)*cos(alpha_13)*sin(alpha_5)*sin(c_) - 2*R*r*cos(a_)*cos(c_)*sin(alpha_5)*sin(alpha_13) + 2*a56*r*cos(alpha_13)*cos(b_)*cos(c_)*sin(alpha_5) + 2*R*r*cos(alpha_5)*cos(b_)*sin(alpha_13)*sin(c_) - 2*r*zc*cos(a_)*cos(alpha_13)*cos(c_)*sin(b_) - 2*a56*r*cos(b_)*sin(alpha_5)*sin(alpha_13)*sin(c_) + 2*r*yc*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) + 2*r*zc*cos(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 2*r*yc*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 2*a56*r*cos(alpha_5)*cos(al\\\r\npha_13)*cos(c_)*sin(a_)*sin(b_) - 2*R*r*cos(alpha_13)*cos(c_)*sin(a_)*sin(alpha_5)*sin(b_) + 2*a56*r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) + 2*R*r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)*sin(c_))^(1/2)
% q2 = (2^(1/2)*(2*R^2 + 2*a56^2 + 2*r^2 + 2*xc^2 + 2*yc^2 + 2*zc^2 - 2*a56*xc*sin(alpha_5) + 2*R*xc*cos(alpha_5) + 2*a56*yc*cos(alpha_5) + 2*R*yc*sin(alpha_5) - 2*3^(1/2)*R*yc*cos(alpha_5) + 2*3^(1/2)*a56*xc*cos(alpha_5) + 2*3^(1/2)*R*xc*sin(alpha_5) + 2*3^(1/2)*a56*yc*sin(alpha_5) - 2*r*xc*cos(alpha_13)*cos(b_)*cos(c_) - 2*r*yc*cos(a_)*cos(alpha_13)*sin(c_) - 2*r*yc*cos(a_)*cos(c_)*sin(alpha_13) + 2*r*xc*cos(b_)*sin(alpha_13)*sin(c_) - 2*r*zc*cos(alpha_13)*sin(a_)*sin(c_) - 2*r*zc*cos(c_)*sin(a_)*sin(alpha_13) - 2*3^(1/2)*r*xc*cos(alpha_13)*cos(b_)*sin(c_) - 2*3^(1/2)*r*xc*cos(b_)*cos(c_)*sin(alpha_13) + 2*3^(1/2)*r*zc*cos(alpha_13)*cos(c_)*sin(a_) - 2*3^(1/2)*r*yc*cos(a_)*sin(alpha_13)*sin(c_) - 2*3^(1/2)*r*zc*sin(a_)*sin(alpha_13)*sin(c_) - 3*R*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_) - R*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*cos(c_) - a56*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*sin(c_) - a56*r*cos(a_)*cos(alpha_5)*cos(c_)*sin(alpha_13) + 3*a56*r*cos(a_)*cos(alpha_13)*cos(c_)*sin(alpha_5\\\r\n) + 3*R*r*cos(a_)*cos(alpha_5)*sin(alpha_13)*sin(c_) - R*r*cos(a_)*cos(alpha_13)*sin(alpha_5)*sin(c_) - R*r*cos(a_)*cos(c_)*sin(alpha_5)*sin(alpha_13) - 3*a56*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(c_) - 3*a56*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_13) + a56*r*cos(alpha_13)*cos(b_)*cos(c_)*sin(alpha_5) + R*r*cos(alpha_5)*cos(b_)*sin(alpha_13)*sin(c_) - 3*R*r*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) - 3*R*r*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13) + 2*r*zc*cos(a_)*cos(alpha_13)*cos(c_)*sin(b_) - 3*a56*r*cos(a_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - a56*r*cos(b_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - 2*r*yc*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) - 2*r*zc*cos(a_)*sin(alpha_13)*sin(b_)*sin(c_) + 2*r*yc*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) + 2*3^(1/2)*r*yc*cos(a_)*cos(alpha_13)*cos(c_) - a56*r*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) + 3*R*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 3*R*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - R*r*cos(alpha_13)*cos(\\\r\nc_)*sin(a_)*sin(alpha_5)*sin(b_) + a56*r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 3*a56*r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) - 3*a56*r*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_) + R*r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)*sin(c_) + 3^(1/2)*a56*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_) + 3^(1/2)*R*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*sin(c_) + 3^(1/2)*R*r*cos(a_)*cos(alpha_5)*cos(c_)*sin(alpha_13) + 3^(1/2)*R*r*cos(a_)*cos(alpha_13)*cos(c_)*sin(alpha_5) - 3^(1/2)*a56*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*cos(c_) - 3^(1/2)*R*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(c_) - 3^(1/2)*R*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_13) - 3^(1/2)*R*r*cos(alpha_13)*cos(b_)*cos(c_)*sin(alpha_5) - 3^(1/2)*a56*r*cos(a_)*cos(alpha_5)*sin(alpha_13)*sin(c_) - 3^(1/2)*a56*r*cos(a_)*cos(alpha_13)*sin(alpha_5)*sin(c_) - 3^(1/2)*a56*r*cos(a_)*cos(c_)*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*R*r*cos(a_)*sin(alpha_5)*sin(alpha_13)*sin(c_) + 3^(1/2)*a56*r*cos(alpha_5)*cos(b\\\r\n_)*sin(alpha_13)*sin(c_) + 3^(1/2)*a56*r*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) + 3^(1/2)*a56*r*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*R*r*cos(b_)*sin(alpha_5)*sin(alpha_13)*sin(c_) + 2*3^(1/2)*r*zc*cos(a_)*cos(alpha_13)*sin(b_)*sin(c_) + 2*3^(1/2)*r*zc*cos(a_)*cos(c_)*sin(alpha_13)*sin(b_) - 2*3^(1/2)*r*yc*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) - 2*3^(1/2)*r*yc*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) + 3^(1/2)*R*r*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) - 3^(1/2)*a56*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) - 3^(1/2)*a56*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - 3^(1/2)*a56*r*cos(alpha_13)*cos(c_)*sin(a_)*sin(alpha_5)*sin(b_) - 3^(1/2)*R*r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 3^(1/2)*R*r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) - 3^(1/2)*R*r*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_) + 3^(1/2)*a56*r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)*sin(c_))^(1/2))/2
% q3 = (2^(1/2)*(2*R^2 + 2*a56^2 + 2*r^2 + 2*xc^2 + 2*yc^2 + 2*zc^2 - 2*a56*xc*sin(alpha_5) + 2*R*xc*cos(alpha_5) + 2*a56*yc*cos(alpha_5) + 2*R*yc*sin(alpha_5) + 2*3^(1/2)*R*yc*cos(alpha_5) - 2*3^(1/2)*a56*xc*cos(alpha_5) - 2*3^(1/2)*R*xc*sin(alpha_5) - 2*3^(1/2)*a56*yc*sin(alpha_5) - 2*r*xc*cos(alpha_13)*cos(b_)*cos(c_) - 2*r*yc*cos(a_)*cos(alpha_13)*sin(c_) - 2*r*yc*cos(a_)*cos(c_)*sin(alpha_13) + 2*r*xc*cos(b_)*sin(alpha_13)*sin(c_) - 2*r*zc*cos(alpha_13)*sin(a_)*sin(c_) - 2*r*zc*cos(c_)*sin(a_)*sin(alpha_13) + 2*3^(1/2)*r*xc*cos(alpha_13)*cos(b_)*sin(c_) + 2*3^(1/2)*r*xc*cos(b_)*cos(c_)*sin(alpha_13) - 2*3^(1/2)*r*zc*cos(alpha_13)*cos(c_)*sin(a_) + 2*3^(1/2)*r*yc*cos(a_)*sin(alpha_13)*sin(c_) + 2*3^(1/2)*r*zc*sin(a_)*sin(alpha_13)*sin(c_) - 3*R*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_) - R*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*cos(c_) - a56*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*sin(c_) - a56*r*cos(a_)*cos(alpha_5)*cos(c_)*sin(alpha_13) + 3*a56*r*cos(a_)*cos(alpha_13)*cos(c_)*sin(alpha_5\\\r\n) + 3*R*r*cos(a_)*cos(alpha_5)*sin(alpha_13)*sin(c_) - R*r*cos(a_)*cos(alpha_13)*sin(alpha_5)*sin(c_) - R*r*cos(a_)*cos(c_)*sin(alpha_5)*sin(alpha_13) - 3*a56*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(c_) - 3*a56*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_13) + a56*r*cos(alpha_13)*cos(b_)*cos(c_)*sin(alpha_5) + R*r*cos(alpha_5)*cos(b_)*sin(alpha_13)*sin(c_) - 3*R*r*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) - 3*R*r*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13) + 2*r*zc*cos(a_)*cos(alpha_13)*cos(c_)*sin(b_) - 3*a56*r*cos(a_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - a56*r*cos(b_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - 2*r*yc*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) - 2*r*zc*cos(a_)*sin(alpha_13)*sin(b_)*sin(c_) + 2*r*yc*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 2*3^(1/2)*r*yc*cos(a_)*cos(alpha_13)*cos(c_) - a56*r*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) + 3*R*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 3*R*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - R*r*cos(alpha_13)*cos(\\\r\nc_)*sin(a_)*sin(alpha_5)*sin(b_) + a56*r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) - 3*a56*r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) - 3*a56*r*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_) + R*r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)*sin(c_) - 3^(1/2)*a56*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*cos(c_) - 3^(1/2)*R*r*cos(a_)*cos(alpha_5)*cos(alpha_13)*sin(c_) - 3^(1/2)*R*r*cos(a_)*cos(alpha_5)*cos(c_)*sin(alpha_13) - 3^(1/2)*R*r*cos(a_)*cos(alpha_13)*cos(c_)*sin(alpha_5) + 3^(1/2)*a56*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*cos(c_) + 3^(1/2)*R*r*cos(alpha_5)*cos(alpha_13)*cos(b_)*sin(c_) + 3^(1/2)*R*r*cos(alpha_5)*cos(b_)*cos(c_)*sin(alpha_13) + 3^(1/2)*R*r*cos(alpha_13)*cos(b_)*cos(c_)*sin(alpha_5) + 3^(1/2)*a56*r*cos(a_)*cos(alpha_5)*sin(alpha_13)*sin(c_) + 3^(1/2)*a56*r*cos(a_)*cos(alpha_13)*sin(alpha_5)*sin(c_) + 3^(1/2)*a56*r*cos(a_)*cos(c_)*sin(alpha_5)*sin(alpha_13) + 3^(1/2)*R*r*cos(a_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - 3^(1/2)*a56*r*cos(alpha_5)*cos(b\\\r\n_)*sin(alpha_13)*sin(c_) - 3^(1/2)*a56*r*cos(alpha_13)*cos(b_)*sin(alpha_5)*sin(c_) - 3^(1/2)*a56*r*cos(b_)*cos(c_)*sin(alpha_5)*sin(alpha_13) - 3^(1/2)*R*r*cos(b_)*sin(alpha_5)*sin(alpha_13)*sin(c_) - 2*3^(1/2)*r*zc*cos(a_)*cos(alpha_13)*sin(b_)*sin(c_) - 2*3^(1/2)*r*zc*cos(a_)*cos(c_)*sin(alpha_13)*sin(b_) + 2*3^(1/2)*r*yc*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 2*3^(1/2)*r*yc*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) - 3^(1/2)*R*r*cos(alpha_5)*cos(alpha_13)*cos(c_)*sin(a_)*sin(b_) + 3^(1/2)*a56*r*cos(alpha_5)*cos(alpha_13)*sin(a_)*sin(b_)*sin(c_) + 3^(1/2)*a56*r*cos(alpha_5)*cos(c_)*sin(a_)*sin(alpha_13)*sin(b_) + 3^(1/2)*a56*r*cos(alpha_13)*cos(c_)*sin(a_)*sin(alpha_5)*sin(b_) + 3^(1/2)*R*r*cos(alpha_5)*sin(a_)*sin(alpha_13)*sin(b_)*sin(c_) + 3^(1/2)*R*r*cos(alpha_13)*sin(a_)*sin(alpha_5)*sin(b_)*sin(c_) + 3^(1/2)*R*r*cos(c_)*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_) - 3^(1/2)*a56*r*sin(a_)*sin(alpha_5)*sin(alpha_13)*sin(b_)*sin(c_))^(1/2))/2

%%
theta_0 = pi/4;
l_0 = 0.1305;
max_iter = 10;
tol = 1e-20;
par_sel = [4,5,6];
ser_sel = [10,11,7];
qp_0 = [theta_0; theta_0; theta_0; l_0; l_0; l_0; l_0*sin(theta_0); 0; 0; 0; 0; 0 ];

% q_ser = [-0.455173816243597, -0.469953566203639, 0.078685071468735 ]';
alpha = -10:10;
beta = -10:10;
x = 0.08:0.001:0.12;

for a = alpha
    for b = beta
        for xc = x
            q_ser = [deg2rad(a), deg2rad(b), xc]
            qp_ik = rps_solve(qp_0,ser_sel,q_ser,max_iter,tol);
            qp_ik_2 = MEII_IK(q_ser);
            abs(qp_ik-qp_ik_2')
        end
    end
end

