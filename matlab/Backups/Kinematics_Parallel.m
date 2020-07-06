clear all; close all; clc;

% Variables that can be controlled
syms q1 q2 q31 q32 q41 q42 q51 q52...
     q1d q2d q31d q32d q41d q42d q51d q52d

qs = 0.11769644556;
 
% Our input variables
variables       = [ q31,  q32,  q41,  q42,  q51,  q52];
variable_values = [-0.453979582395,  qs, -0.453979582395,  qs, -0.453979582395,  qs];

variables_d     = [ q31d,  q32d,  q41d,  q42d,  q51d,  q52d];

R       = 0.1044956;   % meters
alpha_5 = 0.094516665; % radians
a4      = 0.159385;    % meters
a5      = 0.0268986;   % meters
a6      = 0.027282;    % meters
a56     = -(a5-a6);    % meters
alpha_13 = 5*pi/180;

% Number of links in robot
num_of_links = 6;

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
        %3      4     5       6     7       8
a     = [2,     3,    2,      5,    2,      7   ]-2;
mu    = [0,     1,    0,      1,    0,      1   ];
sigma = [0,     1,    0,      1,    0,      1   ];
gamma = [0,     0,    2*pi/3, 0,    -2*pi/3, 0   ];
b     = [0,     0,    0,      0,    0,      0   ];
alpha = [-pi/2, pi/2, -pi/2,  pi/2, -pi/2,  pi/2];
d     = [R,     0,    R,      0,    R,      0   ];
theta = [q31,   0,    q41,    0,    q51,    0   ];
r     = [a56,   q32,  a56,    q42,  a56,    q52 ];

for i = 1:num_of_links
    T(i,:,:) = mdhtf(gamma(i), b(i), alpha(i), d(i), theta(i), r(i));
end

% To get to the position and orientation at each coordinate frame
for i = 1:num_of_links
    if a(i) == 0
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
for i = 1:num_of_links
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


% zp = cross( O(5,:) - O(3,:) , O(7,:) - O(3,:)); % b
% zp = simplify(zp/norm(zp));
% 
% yp = O(7,:)-O(5,:); % (a)
% yp = -simplify(yp/norm(yp));
% % 
% % a = yp;
% % b = zp;
% % 
% % a_par = (a(1)*a(1)+a(2)*a(2)+a(3)*a(3))/(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))*b;
% % a_per = a - a_par;
% % 
% % w = cross(b,a_per);
% % 
% % x1 = cos(-alpha_13)/norm(a_per);
% % x2 = sin(-alpha_13)/norm(w);
% % a_per_rot = norm(a_per)*(x1*a_per+x2*w);
% % a_rot = a_per_rot + a_par;
% % yp = -a_rot/norm(a_rot);
% 
% xp = cross(yp,zp);
% r = 0.05288174521;
% % P = O(3,:)+r*sin(alpha_13)*-yp + r*cos(alpha_13)*-xp;
% % P = O(3,:)+ r*-xp;
% mid2 = (O(5,:) + O(7,:))/2;
% line2 = (O(3,:) - mid2)/norm((O(3,:) - mid2));
% P = O(3,:)-line2*r;

zp = cross( O(5,:) - O(3,:) , O(7,:) - O(3,:)); % b
zp = simplify(zp/norm(zp));
% xp = cross(yp,zp);
r = 0.05288174521;
% P = O(3,:)+r*sin(alpha_13)*-yp + r*cos(alpha_13)*-xp;
% P = O(3,:)+ r*-xp;
mid2 = (O(5,:) + O(7,:))/2;
line2 = (O(3,:) - mid2)/norm((O(3,:) - mid2));
P = O(3,:)-line2*r;

p2 = simplify(expand(rotate_point_vector(O(5,:),zp, P,-alpha_13)));
p3 = simplify(expand(rotate_point_vector(O(7,:),zp, P,-alpha_13)));

yp = (p2-p3)/norm(p2-p3);
xp = cross(yp,zp);


P_eval = subs(P,variables,variable_values);
scatter3(P_eval(1),P_eval(2),P_eval(3));
plot3([P_eval(1), P_eval(1)+0.01*xp(1)],[P_eval(2), P_eval(2)+0.01*xp(2)],[P_eval(3), P_eval(3)+0.01*xp(3)],'r')
plot3([P_eval(1), P_eval(1)+0.01*yp(1)],[P_eval(2), P_eval(2)+0.01*yp(2)],[P_eval(3), P_eval(3)+0.01*yp(3)],'g')
plot3([P_eval(1), P_eval(1)+0.01*zp(1)],[P_eval(2), P_eval(2)+0.01*zp(2)],[P_eval(3), P_eval(3)+0.01*zp(3)],'b')

%% Parallel kinematics
syms px py pz w_alpha w_beta w_gamma r...
     q31 q32 q41 q42 q51 q52

alpha_13 = 5*pi/180;
r = 0.05288174521;
R       = 0.1044956;   % meters
alpha_5 = 0.094516665; % radians
a4      = 0.159385;    % meters
a5      = 0.0268986;   % meters
a6      = 0.027282;    % meters
a56     = -(a5-a6);    % meters

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
 

% vars = [px, py, pz, w_alpha, w_beta, w_gamma];
% var_vals = [0 0 0.1 0 0 0];
% vpa(subs(solve(h(1)==0,q32),vars,var_vals),3)
% vpa(subs(solve(h(2)==0,q42),vars,var_vals),3)
% vpa(subs(solve(h(3)==0,q52),vars,var_vals),3)




% solve([h(1)==0,h(2)==0,h(3)==0],[py,pz,w_gamma])
% py_h1 = vpa(simplify(expand(solve(h(1)==0,py))),2)
% pz_h2 = vpa(simplify(expand(solve(subs(h(2),py,py_h1(1))==0,pz))),2)
% w_gamma = vpa(simplify(expand(solve(subs(h(3),[py,pz],[py_h1(1),pz_h2(1)])==0,w_gamma))),2)

%% Parallel velocity analysis
clear all; close all; clc;

% Variables that can be controlled
syms q1 q2 q31 q32 q41 q42 q51 q52...
     q1d q2d q31d q32d q41d q42d q51d q52d

qs = 0.11769644556;
 
% Our input variables
variables       = [ q31,  q32,  q41,  q42,  q51,  q52];
variable_values = [-0.453979582395,  qs, -0.453979582395,  qs, -0.453979582395,  qs];

R       = 0.1044956;   % meters
alpha_5 = 0.094516665; % radians
a4      = 0.159385;    % meters
a5      = 0.0268986;   % meters
a6      = 0.027282;    % meters
a56     = -(a5-a6);    % meters
alpha_13 = 5*pi/180;

% Number of links in robot
num_of_links = 6;

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
        %3      4     5       6     7       8
a     = [2,     3,    2,      5,    2,      7   ]-2;
mu    = [0,     1,    0,      1,    0,      1   ];
sigma = [0,     1,    0,      1,    0,      1   ];
gamma = [0,     0,    2*pi/3, 0,    4*pi/3, 0   ];
b     = [0,     0,    0,      0,    0,      0   ];
alpha = [-pi/2, pi/2, -pi/2,  pi/2, -pi/2,  pi/2];
d     = [R,     0,    R,      0,    R,      0   ];
theta = [q31,   0,    q41,    0,    q51,    0   ];
r     = [a56,   q32,  a56,    q42,  a56,    q52 ];

for i = 1:num_of_links
    T(i,:,:) = mdhtf(gamma(i), b(i), alpha(i), d(i), theta(i), r(i));
end

% To get to the position and orientation at each coordinate frame
for i = 1:num_of_links
    if a(i) == 0
        Tcom(i,:,:) = T(i,:,:);
    else
        Tcom(i,:,:) = squeeze(Tcom(a(i),:,:))*squeeze(T(i,:,:));
    end
    O(i+1,1) = Tcom(i,1,4);
    O(i+1,2) = Tcom(i,2,4);
    O(i+1,3) = Tcom(i,3,4);
    X(i+1,:) = O(i+1,:)+(squeeze(Tcom(i,1:3,1:3))*[1;0;0])'; 
    Y(i+1,:) = O(i+1,:)+(squeeze(Tcom(i,1:3,1:3))*[0;1;0])';
    Z(i+1,:) = O(i+1,:)+(squeeze(Tcom(i,1:3,1:3))*[0;0;1])';
end

% for i = 1:6
%     if a(i) == 0
%         r_hat = O(i+1,:);
%     else
%         r_hat = O(i+1,:)-O(a(i)+1,:);
%     end
%     r_hat_mat = [0,       -r_hat(3), r_hat(2);
%                     r_hat(3), 0,       -r_hat(1);
%                    -r_hat(2), r_hat(1), 0];
%     R = squeeze(T(i,1:3,1:3)).';
%     j_T_bar_i = [R,         -R*r_hat_mat;
%                  zeros(3,3), R];
% end
% 
% for i = 1:2
%     a_bar = (squeeze(Tcom(i,1:3,1:3))*[0;0;1]);
%     tk{i} = [sigma(i)*a_bar + (1-sigma(i))*(cross(a_bar,(O(3,:)-O(i+1,:)).'));
%              (1-sigma(i))*a_bar];
% end
% 
% recip = [sin(q31) 0 cos(q31) 0 0 0];



% yp = O(7,:)-O(5,:); % (a)
% yp = simplify(yp/norm(yp));
% 
% a = yp;
% b = zp;
% 
% a_par = (a(1)*a(1)+a(2)*a(2)+a(3)*a(3))/(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))*b;
% a_per = a - a_par;
% 
% w = cross(b,a_per);
% 
% x1 = cos(-alpha_13)/norm(a_per);
% x2 = sin(-alpha_13)/norm(w);
% a_per_rot = norm(a_per)*(x1*a_per+x2*w);
% a_rot = a_per_rot + a_par;
% yp = -a_rot/norm(a_rot);

% P = center_triangle(O(3,:), O(5,:), O(7,:));
zp = cross( O(5,:) - O(3,:) , O(7,:) - O(3,:)); % b
zp = simplify(zp/norm(zp));
% xp = cross(yp,zp);
r = 0.05288174521;
% P = O(3,:)+r*sin(alpha_13)*-yp + r*cos(alpha_13)*-xp;
% P = O(3,:)+ r*-xp;
mid2 = (O(5,:) + O(7,:))/2;
line2 = (O(3,:) - mid2)/norm((O(3,:) - mid2));
P = O(3,:)-line2*r;

p2 = simplify(expand(rotate_point_vector(O(5,:),zp, P,alpha_13)));
p3 = simplify(expand(rotate_point_vector(O(7,:),zp, P,alpha_13)));

yp = (p3-p2)/norm(p3-p2);
xp = cross(yp,zp);

P_eval = subs(P,variables,variable_values);

%%
clc; clear all;
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z r alpha_13
O1 = [p1x, p1y, p1z]; O2 = [p2x, p2y, p2z]; O3 = [p3x, p3y, p3z];

zp = cross( O2 - O1 , O3 - O1); % b
zp = simplify(zp/norm(zp));

mid2 = (O2 + O3)/2;
line2 = (O1 - mid2)/norm((O1 - mid2));
P = O1-line2*r;

p2 = rotate_point_vector(O2,zp, P,-alpha_13);
p3 = rotate_point_vector(O3,zp, P,-alpha_13);

yp = (p2-p3)/norm(p2-p3);
xp = cross(yp,zp);

R = [xp.', yp.', zp.'];
% P_eval = subs(P,variables,variable_values);

%%
syms px py pz alpha beta gamma

zp = cross( O2 - O1 , O3 - O1); % b
zp = simplify(zp/norm(zp));

mid2 = (O2 + O3)/2;
line2 = (O1 - mid2)/norm((O1 - mid2));
P = [px, py, pz];

p2 = rotate_point_vector(O2,zp, P,-alpha_13);
p3 = rotate_point_vector(O3,zp, P,-alpha_13);

yp = (p2-p3)/norm(p2-p3);
xp = cross(yp,zp);

R = [xp.', yp.', zp.'];

Rx(alpha)*Ry(beta)*Rz(gamma)