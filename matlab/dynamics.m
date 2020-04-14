clc; clear all; close all;

% Variables that can be controlled
syms q11 q21 q12 q22 q13 q23...
     q11d q21d q12d q22d q13d q23d...
     vpx vpy vpz wpx wpy wpz...
     vpx_dot vpy_dot vpz_dot wpx_dot wpy_dot wpz_dot...
     px  py  pz  a b c...
     d1 dp
     
gamma = [0, 2*pi/3, 4*pi/3;
         0, 0,      0];
 
q = [q11 q12 q13; 
     q21 q22 q23];

q_dot = [q11d q12d q13d; 
         q21d q22d q23d];
     
qa = [q11, q12, q13].';

qa_dot = [q11d, q12d, q13d].';

%%%%%%%%% KINEMATIC MODEL %%%%%%%%%%%%%

R_0_P = Rx(a)*Ry(b)*Rz(c);
P_0 = [px, py, pz].';

% Eqn 5 NOT SURE ABOUT THIS
T_0_P = [R_0_P(1:3,1:3), P_0;
         0, 0, 0, 1];

for i = 1:3
    % Eqn 8
    J_1i_{i} = [-q(2,i), 0;
                0,      1;
                0,      0];
    
    % Eqn 9
    J_1i_di{i} = [-q(2,i), 0;
                   0,      1];
            
    % Eqn 11
    J_1i_di_inv{i} = [-1/q(2,i), 0;
                       0,        1];
            
    % Eqn 12 (intermediate calc)
    R_1_0{i} = [ cos(gamma(1,i))*cos(q(1,i)),  sin(gamma(1,i))*cos(q(1,i)), -sin(q(1,i));
                -cos(gamma(1,i))*sin(q(1,i)), -sin(gamma(1,i))*sin(q(1,i)), -cos(q(1,i));
                -sin(gamma(1,i)),              cos(gamma(1,i)),               0        ];
    
    % THIS CHANGES FOR OUR MECHANISM -> 0 becomes a56
    % Eqn 3
    B_0_{i} = [cos(gamma(1,i))*d1, sin(gamma(1,i))*d1, 0].';
    
    % Eqn 4
    P_0_{i} = [cos(gamma(1,i))*dp, sin(gamma(1,i))*dp, 0].';
    
    % Eqn 12 (intermediate calc)
    PP_0_{i} = P_0_{i} - P_0;
 
    % Eqn 12 (intermediate calc)
    PP_hat_0_{i} = crossprodmat(PP_0_{i});
    
    % Eqn 12 (intermediate calc)
    J_rhs = R_1_0{i}*PP_hat_0_{i};
    
    % Eqn 12
    J_1i_v{i} = [R_1_0{i}(1,:), -J_rhs(1,:);
                R_1_0{i}(2,:), -J_rhs(2,:)];
end

% Eqn 17
J_0_p_inv = [J_1i_di_inv{1}(2,:)*J_1i_v{1};
             J_1i_di_inv{2}(2,:)*J_1i_v{2};
             J_1i_di_inv{3}(2,:)*J_1i_v{3}];

v_0_p = [vpx vpy vpz].';
w_0_p = [wpx wpy wpz].';

% V_0_p = [v_0_p.' w_0_p.'].';

% Eqn 21 THIS IS WRONG
B = [R_1_0{1}(3,:);
     R_1_0{2}(3,:);
     R_1_0{3}(3,:)];

% Eqn 20 (intermediate calc)
A_temp{1} = R_1_0{1}*PP_hat_0_{1};
A_temp{2} = R_1_0{2}*PP_hat_0_{2};
A_temp{3} = R_1_0{3}*PP_hat_0_{3};

% Eqn 20
A = [A_temp{1}(3,:);
     A_temp{2}(3,:);
     A_temp{3}(3,:)];

% Eqn 23
C = A\B;

% Eqn 25
a_0_r = [eye(3);
         C];

% Eqn 26
J_0_r_inv = J_0_p_inv*a_0_r;

% Eqn 27 (intermediate calc)
J_0_r = inv(J_0_r_inv);

% Eqn 27
% v_0_p = J_0_r*qa_dot;

% Eqn 28
% w_0_p = C*J_0_r*qa_dot;

% Eqn 29
J_0_p = [J_0_r;
         C*J_0_r];

%% Dynamic Model of the Links

[M,V,G] = LinkDynamics();

%% Inverse Dynamic Model

syms Mp MXp MYp MZp...
Ipxx Ipyy Ipzz Ipxy Ipxz Ipyz...
H11 H12 H13 H21 H22 H23...
T1 T2 T3...
g_

g = [g_, 0, 0].';

MSp = [MXp MYp MZp].';
Ip = [Ipxx -Ipxy -Ipxz;
     -Ipxy Ipyy -Ipyz;
     -Ipxz -Ipyz Ipxx];
 
vp = [vpx vpy vpz].';
wp = [wpx wpy wpz].';
 
% vp_dot = [vpx_dot vpy_dot vpz_dot].';
% wp_dot = [wpx_dot wpy_dot wpz_dot].';

Jp = [Mp*eye(3),        -crossprodmat(MSp);
      crossprodmat(MSp) Ip];
  
for i = 1:3
    % Eqn 32
    PP_0_dot_{i} = cross(w_0_p,PP_0_{i}); 
end

% Eqn 36
% REDUCE THIS TO QUANTITIES - mainly the sqrt(3)/2s and 1/2s
A_dot = [PP_0_dot_{1}(3),             0,                   -PP_0_dot_{1}(1);
         -sqrt(3)/2*PP_0_dot_{2}(3), -1/2*PP_0_dot_{2}(3), 1/2*PP_0_dot_{2}(2)+sqrt(3)/2*PP_0_dot_{2}(1);
         sqrt(3)/2*PP_0_dot_{3}(3),  -1/2*PP_0_dot_{3}(3), 1/2*PP_0_dot_{3}(2)-sqrt(3)/2*PP_0_dot_{3}(1)];
  
% Eqn 39
a_dot_0_r = [zeros(3,3);-inv(A)*A_dot*C];
  
% Fp = Jp*[(vp_dot - g).', wp_dot.'].'+[cross(wp,cross(wp,MSp)).', cross(wp,Ip*wp).'].';

sum_JAJ = zeros(6,6);
sum_JAJavJqh = zeros(6,1);

for i = 1:3
    
    J_1i_dot{i} = [-q_dot(2,i), 0;
                    0,          1;
                    0,          0];
                
    R_1_0_dot{i} = [-sin(q(1,i))*cos(gamma(1,i))*q_dot(1,i), -sin(q(1,i))*sin(gamma(1,i))*q_dot(1,i), -cos(q(1,i))*q_dot(1,i);
                    -cos(q(1,i))*cos(gamma(1,i))*q_dot(1,i), -cos(q(1,i))*sin(gamma(1,i))*q_dot(1,i),  sin(q(1,i))*q_dot(1,i);
                     0,                                       0,                                       0];
                 
    PP_hat_0_dot{i} = [  0,  vpz, -vpy;
                      -vpz,    0,  vpx;
                       vpy, -vpx,    0];

    RPP_dot{i} = -(R_1_0{i}*PP_hat_0_dot{i} + R_1_0_dot{i}*PP_hat_0_{i});
                   
    J_1i_vi_dot{i} = [R_1_0_dot{i}(1,:), RPP_dot{i}(1,:);
                      R_1_0_dot{i}(2,:), RPP_dot{i}(2,:)];
    
    Ai{i} = M(1+2*(i-1):2+2*(i-1),1+2*(i-1):2+2*(i-1));
    hi{i} = V(1+2*(i-1):2+2*(i-1))+G(1+2*(i-1):2+2*(i-1));
    
    Axi{i} = J_1i_di_inv{i}.'*Ai{i}*J_1i_di_inv{i};
    hxi{i} = J_1i_di_inv{i}.'*hi{i};
    
    sum_JAJ = sum_JAJ + J_1i_v{i}.'*Axi{i}*J_1i_v{i};
    sum_JAJavJqh = sum_JAJavJqh + J_1i_v{i}.'*(Axi{i}*(J_1i_vi_dot{i}*a_0_r*vp-J_1i_dot{i}(1:2,:)*q_dot(1:2,i))+hxi{i});
end

%NOT SURE ABOUT JP VS IP HERE AND IN NEXT EQUATION
A_robot = a_0_r.'*(Jp + sum_JAJ)*a_0_r;

h_robot = a_0_r.'*(Jp + sum_JAJ)*a_dot_0_r*vp + ...
          a_0_r.'*([cross(wp,cross(wp,MSp)); cross(wp,Ip*wp)]-[Mp*eye(3);crossprodmat(MSp)]*g)+ ...
          a_0_r.'*sum_JAJavJqh;

Gamma = [T1, T2, T3].';

vp_dot = inv(A_robot)*(J_0_r_inv.'*Gamma - h_robot);

%% %%%%%%%%% SECOND ORDER KINEMATIC MODEL %%%%%%%%%%%%%
     
% for i = 1:3
%     % Eqn 32
%     J_1i_di_dot{i} = [-q_dot(2,i)*q_dot(1,i),0].';
% end

% Eqn 34
V_dot_0_p = a_0_r*vp_dot + a_dot_0_r*v_0_p;

simplified = simplify(expand(V_dot_0_p(5)))

% for i = 1:3
%     % Eqn 33
%     v_dot_1i_{i} = R_1_0{i}*([eye(3), crossprodmat(PP_hat_0_{i})]*V_dot_0_p + cross(w_0_p,cross(w_0_p,PP_0_{i})));
%     
%     qdd(1:2,i) = J_1i_di_inv{i}*(v_dot_1i_{i}(1:2) - J_1i_di_dot{i});
% end