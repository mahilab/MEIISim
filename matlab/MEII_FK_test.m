% clc; clear all;
% syms t1 t2 t3 l1 l2 l3 p
% 
% % variables = [l1 l2 l3];
% % values = [0.1305 0.1305 0.1305];
% r = 0.05288174521; R = 0.10449630336;
% l1 = 0.1305/R; l2 = 0.1305/R; l3 = 0.1305/R;
% p = r/R;
% theta1 = pi/4; theta2 = pi/4; theta3 = pi/4;
% % t1 = 2.73495967-pi/2; t2 = 2.73495967-pi/2; t3 = 2.73495967-pi/2;
% thetas = [theta1, theta2, theta3].';
% 
% eq(1,1) = l1^2 + l2^2 + 3 - 3*p^2 + l1*l2*cos(t1)*cos(t2) - 2*l1*l2*sin(t1)*sin(t2) - 3*l1*cos(t1)-3*l2*cos(t2);
% eq(2,1) = l2^2 + l3^2 + 3 - 3*p^2 + l2*l3*cos(t2)*cos(t3) - 2*l2*l3*sin(t2)*sin(t3) - 3*l2*cos(t2)-3*l3*cos(t3);
% eq(3,1) = l3^2 + l1^2 + 3 - 3*p^2 + l3*l1*cos(t3)*cos(t1) - 2*l3*l1*sin(t3)*sin(t1) - 3*l3*cos(t3)-3*l1*cos(t1);
% 
% jac = jacobian(eq,[t1,t2,t3]);
% old_thetas = [inf, inf, inf].';
% iter = 0;
% while abs(rms(old_thetas-thetas)) > 1e-50 && iter < 20
%     old_thetas = thetas;
%     thetas = vpa(thetas - subs(jac,[t1,t2,t3].',thetas)\subs(eq,[t1 t2 t3].',thetas),32)
%     iter = iter + 1
% end
% 
% % t2_solved = simplify(expand(solve(eq1==0,t2,'Real',true)))
% % t3_solved = simplify(solve(subs(eq2==0,t2,t2_solved(1)),t3,'Real',true))
% % t1_solved = solve(subs(eq3==0,t3,t3_solved(1)),t1)
% 
% % vpa(subs(t1_solved,variables,values),12)
% vpa(2.73495967-pi/2-thetas(1),12)

%%
syms n1 n2 n3 o1 o2 o3 a1 a2 a3 R r a56 alpha_5 alpha_13 l1 l2 l3 theta1 theta2 theta3
  
T = [n1 o1 a1;
     n2 o2 a2;
     n3 o3 a3];

R = 0.1044956; % [m]
r = 0.05288174521; % [m]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]
alpha_5 = 0.094516665; % [rad]
alpha_13 = 5*pi/180; % [rad]

gamma = [0, -2*pi/3, 2*pi/3];

theta = [1.06554686464; 1.0871646663; 0.582798766322];
l = [0.119273272676; 0.101143787974; 0.0837368660614];
% theta = ones(3,1)*pi-(1.164174936153484886976163989065+pi/2);

% l = [l1 l2 l3];
% theta = [theta1 theta2 theta3];

for i = 1:3
    B_temp = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-theta(i))*TRANSx(l(i));
    B{i} = B_temp(1:3,4);
    b_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    b{i} = b_temp{i}(1:3,4);
    
end

P_c = (B{1} + B{2} + B{3})/3
P_c = [px py pz].';
for i = 1:2
    eqs_(-2+3*i:0+3*i) = T*b{i} + P_c == B{i};
end

[A_lin,B_lin] = equationsToMatrix(eqs_, [a1 a2 a3 o1 o2 o3]);

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
% n3_ =  o1_*a2_ - o2_*a1_;

beta = vpa(asin(n2_),4)
alpha = vpa(acos(n1_/cos(beta)),4)
gamma = vpa(acos(o2_/cos(beta)),4)

%%
% l1 = 0.1305; l2 = 0.1305; l3 = 0.1305;
l1 = 0.119273272676; l2 = 0.101143787974; l3 = 0.0837368660614;
% l1 = l1/R; l2 = l2/R; l3 = l3/R;
R = 0.1044956; % [m]
r = 0.05288174521; % [m]
alpha_5 = 0.094516665; % [rad]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]

thetas = [pi/4,pi/4,pi/4].';
old_thetas = [pi pi pi].';
iter = 0

while abs(rms(old_thetas-thetas)) > 1e-20 && iter < 20
    old_thetas = thetas;
    t1 = thetas(1); t2 = thetas(2); t3 = thetas(3);
    jac = [ -(l1*(3*R*cos(t1) - l2*cos(t1)*sin(t2) - 2*l2*cos(t2)*sin(t1) + 3^(1/2)*a56*cos(t1)))/(2*(abs(l1*cos(t1) - l2*cos(t2))^2 + abs((3*R*cos(alpha_5))/2 + (3*a56*sin(alpha_5))/2 - l1*cos(alpha_5)*sin(t1) - (l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l2*sin(alpha_5)*sin(t2))/2 - (3^(1/2)*R*cos(alpha_5))/2 - (3^(1/2)*a56*sin(alpha_5))/2 + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2)^2)^(1/2)), (l2*(2*l1*cos(t1)*sin(t2) - 3*R*cos(t2) + l1*cos(t2)*sin(t1) + 3^(1/2)*a56*cos(t2)))/(2*(abs(l1*cos(t1) - l2*cos(t2))^2 + abs((3*R*cos(alpha_5))/2 + (3*a56*sin(alpha_5))/2 - l1*cos(alpha_5)*sin(t1) - (l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l2*sin(alpha_5)*sin(t2))/2 - (3^(1/2)*R*cos(alpha_5))/2 - (3^(1/2)*a56*sin(alpha_5))/2 + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2)^2)^(1/2)),0;0,-(l2*(3*R*cos(t2) - l3*cos(t2)*sin(t3) - 2*l3*cos(t3)*sin(t2) + 3^(1/2)*a56*cos(t2)))/(2*(abs(l2*cos(t2) - l3*cos(t3))^2 + abs((l3*cos(alpha_5)*sin(t3))/2 - (l2*cos(alpha_5)*sin(t2))/2 + 3^(1/2)*a56*cos(alpha_5) - 3^(1/2)*R*sin(alpha_5) + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((l2*sin(alpha_5)*sin(t2))/2 - (l3*sin(alpha_5)*sin(t3))/2 - 3^(1/2)*R*cos(alpha_5) - 3^(1/2)*a56*sin(alpha_5) + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2)),  (l3*(2*l2*cos(t2)*sin(t3) - 3*R*cos(t3) + l2*cos(t3)*sin(t2) + 3^(1/2)*a56*cos(t3)))/(2*(abs(l2*cos(t2) - l3*cos(t3))^2 + abs((l3*cos(alpha_5)*sin(t3))/2 - (l2*cos(alpha_5)*sin(t2))/2 + 3^(1/2)*a56*cos(alpha_5) - 3^(1/2)*R*sin(alpha_5) + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((l2*sin(alpha_5)*sin(t2))/2 - (l3*sin(alpha_5)*sin(t3))/2 - 3^(1/2)*R*cos(alpha_5) - 3^(1/2)*a56*sin(alpha_5) + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2));(l1*(l3*cos(t1)*sin(t3) - 3*R*cos(t1) + 2*l3*cos(t3)*sin(t1) + 3^(1/2)*a56*cos(t1)))/(2*(abs(l1*cos(t1) - l3*cos(t3))^2 + abs(l1*cos(alpha_5)*sin(t1) - (3*a56*sin(alpha_5))/2 - (3*R*cos(alpha_5))/2 + (l3*cos(alpha_5)*sin(t3))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l3*sin(alpha_5)*sin(t3))/2 + (3^(1/2)*R*cos(alpha_5))/2 + (3^(1/2)*a56*sin(alpha_5))/2 - (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2)),0, -(l3*(3*R*cos(t3) - 2*l1*cos(t1)*sin(t3) - l1*cos(t3)*sin(t1) + 3^(1/2)*a56*cos(t3)))/(2*(abs(l1*cos(t1) - l3*cos(t3))^2 + abs(l1*cos(alpha_5)*sin(t1) - (3*a56*sin(alpha_5))/2 - (3*R*cos(alpha_5))/2 + (l3*cos(alpha_5)*sin(t3))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l3*sin(alpha_5)*sin(t3))/2 + (3^(1/2)*R*cos(alpha_5))/2 + (3^(1/2)*a56*sin(alpha_5))/2 - (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2))];
    eq =  [(abs(l1*cos(t1) - l2*cos(t2))^2 + abs((3*R*cos(alpha_5))/2 + (3*a56*sin(alpha_5))/2 - l1*cos(alpha_5)*sin(t1) - (l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l2*sin(alpha_5)*sin(t2))/2 - (3^(1/2)*R*cos(alpha_5))/2 - (3^(1/2)*a56*sin(alpha_5))/2 + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2)^2)^(1/2) - 3^(1/2)*r;(abs(l2*cos(t2) - l3*cos(t3))^2 + abs((l3*cos(alpha_5)*sin(t3))/2 - (l2*cos(alpha_5)*sin(t2))/2 + 3^(1/2)*a56*cos(alpha_5) - 3^(1/2)*R*sin(alpha_5) + (3^(1/2)*l2*sin(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((l2*sin(alpha_5)*sin(t2))/2 - (l3*sin(alpha_5)*sin(t3))/2 - 3^(1/2)*R*cos(alpha_5) - 3^(1/2)*a56*sin(alpha_5) + (3^(1/2)*l2*cos(alpha_5)*sin(t2))/2 + (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2) - 3^(1/2)*r; (abs(l1*cos(t1) - l3*cos(t3))^2 + abs(l1*cos(alpha_5)*sin(t1) - (3*a56*sin(alpha_5))/2 - (3*R*cos(alpha_5))/2 + (l3*cos(alpha_5)*sin(t3))/2 + (3^(1/2)*a56*cos(alpha_5))/2 - (3^(1/2)*R*sin(alpha_5))/2 + (3^(1/2)*l3*sin(alpha_5)*sin(t3))/2)^2 + abs((3*a56*cos(alpha_5))/2 - (3*R*sin(alpha_5))/2 + l1*sin(alpha_5)*sin(t1) + (l3*sin(alpha_5)*sin(t3))/2 + (3^(1/2)*R*cos(alpha_5))/2 + (3^(1/2)*a56*sin(alpha_5))/2 - (3^(1/2)*l3*cos(alpha_5)*sin(t3))/2)^2)^(1/2) - 3^(1/2)*r];
    thetas = thetas - jac\eq;
    iter = iter + 1;
end

[pi/2 pi/2 pi/2].' - thetas;
% vpa(pi-2.73495967-thetas(1),12)
%% FSOLVE METHOD TESTING
theta0 = [pi/3,pi/3,pi/3].';
l1 = 0.119273272676; l2 = 0.101143787974; l3 = 0.0837368660614;
l = [l1 l2 l3];
f = @(theta) eqns_solve(theta,l);
[out,fval] = fsolve(f,theta0);

% out
[pi/2 pi/2 pi/2].' - out

%%

clc; clear all;
syms t1 t2 t3 l1 l2 l3 p r R l1 l2 l3 alpha_5 a56

% variables = [l1 l2 l3];
% values = [0.1305 0.1305 0.1305];
% r = 0.05288174521; R = 0.1044956;
% l1 = 0.1305; l2 = 0.1305; l3 = 0.1305;
% alpha_5 = 0.094516665;
% theta1 = pi/4; theta2 = pi/4; theta3 = pi/4;
% a5 = 0.0268986; % [m]
% a6 = 0.0272820; % [m]
% a56 = -(a5-a6); % [m]
% t1 = pi-2.73495967; t2 = pi-2.73495967; t3 = pi-2.73495967;
% thetas = [theta1, theta2, theta3].';
theta = [t1 t2 t3].';
l = [l1 l2 l3].';
gamma = [0, -2*pi/3, 2*pi/3];

% theta = [pi/2-1.0723 pi/2-1.0779 pi/2-0.7026];
% theta = [pi/2-1.1642, pi/2-1.1642, pi/2-1.1642];
% l     = [0.1193 0.1011 0.0837];
% l     = [0.1305 0.1305 0.1305];
% theta   = [pi/2-1.0655 pi/2-1.0872 pi/2-0.5828];
% l       = [0.1193 0.1011 0.0837];

for i = 1:3
    B_temp = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-theta(i))*TRANSx(l(i));
    B{i} = B_temp(1:3,4);
end

% P_c = (B{1}+B{2}+B{3})/3

eq(1,1) = simplify(expand(sqrt((B{1}(1)-B{2}(1))^2+(B{1}(2)-B{2}(2))^2+(B{1}(3)-B{2}(3))^2)-sqrt(3)*r));
eq(2,1) = simplify(expand(sqrt((B{3}(1)-B{2}(1))^2+(B{3}(2)-B{2}(2))^2+(B{3}(3)-B{2}(3))^2)-sqrt(3)*r));
eq(3,1) = simplify(expand(sqrt((B{3}(1)-B{1}(1))^2+(B{3}(2)-B{1}(2))^2+(B{3}(3)-B{1}(3))^2)-sqrt(3)*r));
% eq(1,1) = l1^2 + l2^2 + 3 - 3*p^2 + l1*l2*cos(t1)*cos(t2) - 2*l1*l2*sin(t1)*sin(t2) - 3*l1*cos(t1)-3*l2*cos(t2);
% eq(2,1) = l2^2 + l3^2 + 3 - 3*p^2 + l2*l3*cos(t2)*cos(t3) - 2*l2*l3*sin(t2)*sin(t3) - 3*l2*cos(t2)-3*l3*cos(t3);
% eq(3,1) = l3^2 + l1^2 + 3 - 3*p^2 + l3*l1*cos(t3)*cos(t1) - 2*l3*l1*sin(t3)*sin(t1) - 3*l3*cos(t3)-3*l1*cos(t1);

jac = jacobian(eq,[t1,t2,t3]);
jac = simplify(expand(jac));
old_thetas = [inf, inf, inf].';
iter = 0;
while abs(rms(old_thetas-thetas)) > 1e-50 && iter < 20
    old_thetas = thetas;
    thetas = vpa(thetas - subs(jac,[t1,t2,t3].',thetas)\subs(eq,[t1 t2 t3].',thetas),32)
    iter = iter + 1
end

% t2_solved = simplify(expand(solve(eq1==0,t2,'Real',true)))
% t3_solved = simplify(solve(subs(eq2==0,t2,t2_solved(1)),t3,'Real',true))
% t1_solved = solve(subs(eq3==0,t3,t3_solved(1)),t1)

% vpa(subs(t1_solved,variables,values),12)
vpa(pi-2.73495967-thetas(1),12)

%%
alpha = -25:0.5:25;
beta = -25:0.5:25;
x = 0.08:0.002:0.13;

diffs = zeros(12,size(alpha,2)*size(beta,2)*size(x,2));

i = 0;
incorr = 0;
for a = alpha
    for b = beta
        for xc = x
            q_ik = MEII_IK_3([deg2rad(a),deg2rad(b),x]);
            q_par = q_ik(4:6);
            q_fk = MEII_FK_2(q_par);
            diffs(1:12,i+1) = q_fk-q_ik;
            if abs(sum(diffs(1:12,i+1))) > 1e-5
%                 incorr = incorr + 1;
%                 [a b xc]
%                 [q_ik q_fk diffs(1:12,i+1)]
            end
            i = i + 1;
        end
    end
end

mean_diffs = mean(diffs,2)
max_diffs = max(diffs,[],2)
