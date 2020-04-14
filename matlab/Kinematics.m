clear all; close all; clc;

% Variables that can be controlled
syms q1 q2 q31 q32 q41 q42 q51 q52

% Our input variables
variables       = [q1, q2,  q31,  q32,  q41,  q42,  q51,  q52];
variable_values = [0,  0,  -0.1,  0.1, -0.1,  0.1, -0.1,  0.1];

R       = 0.1044956;   % meters
alpha_5 = 0.094516665; % radians
a4      = 0.159385;    % meters
a5      = 0.0268986;   % meters
a6      = 0.027282;    % meters
a56     = -(a5-a6);    % meters

% Number of links in robot
num_of_links = 8;

% Initialize transformation matrices and position vectors
T    = sym(zeros(num_of_links,4,4));
Tcom = sym(zeros(num_of_links,4,4));
O    = sym(zeros(num_of_links+1,3));
X    = sym(zeros(num_of_links+1,3));
Y    = sym(zeros(num_of_links+1,3));
Z    = sym(zeros(num_of_links+1,3));

% Initialize the first coordinate frame unit vectors
X(1,:) = [0.01,0,0];
Y(1,:) = [0,0.01,0];
Z(1,:) = [0,0,0.01];

% DH Paramaters for each link
        %1   2           3      4     5       6     7       8
a     = [0,  1,          2,     3,    2,      5,    2,      7   ];
mu    = [1,  1,          0,     1,    0,      1,    0,      1   ];
sigma = [0,  0,          0,     1,    0,      1,    0,      1   ];
gamma = [0,  0,          0,     0,    2*pi/3, 0,    4*pi/3, 0   ];
b     = [0,  0,          0,     0,    0,      0,    0,      0   ];
alpha = [0,  pi/2,       -pi/2, pi/2, -pi/2,  pi/2, -pi/2,  pi/2];
d     = [0,  0,          R,     0,    R,      0,    R,      0   ];
theta = [q1, q2+alpha_5, q31,   0,    q41,    0,    q51,    0   ];
r     = [0,  a4,         a56,   q32,  a56,    q42,  a56,    q52 ];

for i = 1:num_of_links
    T(i,:,:) = mdhtf(gamma(i), b(i), alpha(i), d(i), theta(i), r(i));
end

% To get to the position and orientation at each coordinate frame
for i = 1:num_of_links
    if i == 1
        Tcom(i,:,:) = T(i,:,:);
    else
        Tcom(i,:,:) = squeeze(Tcom(a(i),:,:))*squeeze(T(i,:,:));
    end
    O(i+1,1) = subs(Tcom(i,1,4),variables,variable_values);
    O(i+1,2) = subs(Tcom(i,2,4),variables,variable_values);
    O(i+1,3) = subs(Tcom(i,3,4),variables,variable_values);
    X(i+1,:) = O(i+1,:)+subs((squeeze(Tcom(i,1:3,1:3))*[0.01;0;0])',variables,variable_values); 
    Y(i+1,:) = O(i+1,:)+subs((squeeze(Tcom(i,1:3,1:3))*[0;0.01;0])',variables,variable_values);
    Z(i+1,:) = O(i+1,:)+subs((squeeze(Tcom(i,1:3,1:3))*[0;0;0.01])',variables,variable_values);
end

% Plot the origins first with dashed connection lines
for i = 2:num_of_links
    plot3([O(a(i)+1,1),O(i+1,1)],[O(a(i)+1,2),O(i+1,2)],[O(a(i)+1,3),O(i+1,3)],'--ko')
    hold on
end

% Plot each of the coordinate frames
for i = 1:num_of_links+1
    plot3([O(i,1),X(i,1)],[O(i,2),X(i,2)],[O(i,3),X(i,3)],'r')
    hold on
    plot3([O(i,1),Y(i,1)],[O(i,2),Y(i,2)],[O(i,3),Y(i,3)],'g')
    hold on 
    plot3([O(i,1),Z(i,1)],[O(i,2),Z(i,2)],[O(i,3),Z(i,3)],'b')
    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')

% view(gca,90,0)   % (y,z)

% view(gca,-90,90) % (y,x)

view(gca,0,0)    % (-z,x)
camroll(90)

%% Parallel kinematics
syms px py pz w_alpha w_beta w_gamma r

alpha_13 = 5*pi/180;
r = 0.05288174521;

r_p = [px py pz].';

R_5_13 = Rz(-alpha_5)*Ry(w_alpha)*Rz(w_beta)*Rx(w_gamma); % PROBABLY NOT RIGHT
R_5_13 = R_5_13(1:3,1:3);

r_A11 = TRANSx(R)*TRANSy(a56);
r_A11 = r_A11(1:3,4);
r_A21 = Rz(2*pi/3)*TRANSx(R)*TRANSy(a56);
r_A21 = r_A21(1:3,4);
r_A31 = Rz(4*pi/3)*TRANSx(R)*TRANSy(a56);
r_A31 = r_A31(1:3,4);

p_r_A16 = Rz(alpha_13)*TRANSx(r);
p_r_A16 = p_r_A16(1:3,4);
p_r_A26 = Rz(alpha_13+2*pi/3)*TRANSx(r);
p_r_A26 = p_r_A26(1:3,4);
p_r_A36 = Rz(alpha_13+4*pi/3)*TRANSx(r);
p_r_A36 = p_r_A36(1:3,4);

p_r = [px, py, pz].';

r_A16 = p_r + R_5_13.'*p_r_A16;
r_A26 = p_r + R_5_13.'*p_r_A26;
r_A36 = p_r + R_5_13.'*p_r_A36;

h = [(r_A16(1) - r_A11(1))^2 + (r_A16(2) - r_A11(2))^2 + (r_A16(3) - r_A11(3))^2 - q32^2;
     (r_A26(1) - r_A21(1))^2 + (r_A26(2) - r_A21(2))^2 + (r_A26(3) - r_A21(3))^2 - q42^2;
     (r_A36(1) - r_A31(1))^2 + (r_A36(2) - r_A31(2))^2 + (r_A36(3) - r_A31(3))^2 - q52^2];
 
vars = [px, py, pz, w_alpha, w_beta, w_gamma];
var_vals = [0 0 0.1 0 0 0];
vpa(subs(solve(h(1)==0,q32),vars,var_vals),3)
vpa(subs(solve(h(2)==0,q42),vars,var_vals),3)
vpa(subs(solve(h(3)==0,q52),vars,var_vals),3)
% solve(h(1),q32)
% solve(h(1),q32)

%% Parallel velocity analysis

for i = 1:2
    a(i) = 
end