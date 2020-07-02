clear all;

fprintf("Initializing symbolic variables\n");

syms qe qf l1 l2 l3 theta1 theta2 theta3...
     qe_dot qf_dot l1_dot l2_dot l3_dot theta1_dot theta2_dot theta3_dot...
     qe_d_dot qf_d_dot l1_d_dot l2_d_dot l3_d_dot theta1_d_dot theta2_d_dot theta3_d_dot...
     qet(t) qft(t) l1t(t) l2t(t) l3t(t) theta1t(t) theta2t(t) theta3t(t)...
     qe_dott(t) qf_dott(t) l1_dott(t) l2_dott(t) l3_dott(t) theta1_dott(t) theta2_dott(t) theta3_dott(t)...
     wx wy wz vx vy vz px py pz...
     r R alpha_5 a4 a56 alpha_13 l_offset...
     T1 T2 F3 F4 F5...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     Ixx_p Iyy_p Izz_p Mp...
     Ixx_l Iyy_l Izz_l Ml...
     Ixx_f Iyy_f Izz_f Mf...
     Ixx_e Iyy_e Izz_e Me...
     g...
     P_p_x P_p_y P_p_z R_p_x R_p_y R_p_z...
     V_p_x V_p_y V_p_z w_p_x w_p_y w_p_z...
     A_p_x A_p_y A_p_z a_p_x a_p_y a_p_z...
     P_p_xt(t) P_p_yt(t) P_p_zt(t) R_p_xt(t) R_p_yt(t) R_p_zt(t)...
     V_p_xt(t) V_p_yt(t) V_p_zt(t) w_p_xt(t) w_p_yt(t) w_p_zt(t)
 
T = [n1 o1 a1 px;
     n2 o2 a2 py;
     n3 o3 a3 pz;
      0  0  0  1];
 
Tinv = [T(1:3,1:3).',-T(1:3,1:3).'*T(1:3,4);T(4,:)];
 
gamma = [0, -2*pi/3, 2*pi/3];

qs       = [      qe,       qf,       l1,       l2,       l3,       theta1,       theta2,       theta3, P_p_x, P_p_y, P_p_z, R_p_x, R_p_y, R_p_z];
q_dots   = [  qe_dot,   qf_dot,   l1_dot,   l2_dot,   l3_dot,   theta1_dot,   theta2_dot,   theta3_dot, V_p_x, V_p_y, V_p_z, w_p_x, w_p_y, w_p_z];
q_d_dots = [qe_d_dot, qf_d_dot, l1_d_dot, l2_d_dot, l3_d_dot, theta1_d_dot, theta2_d_dot, theta3_d_dot, A_p_x, A_p_y, A_p_z, a_p_x, a_p_y, a_p_z];

qst     = [      qet,      qft,      l1t,      l2t,      l3t,      theta1t,      theta2t,      theta3t, P_p_xt, P_p_yt, P_p_zt, R_p_xt, R_p_yt, R_p_zt];
q_dotst = [  qe_dott,  qf_dott,  l1_dott,  l2_dott,  l3_dott,  theta1_dott,  theta2_dott,  theta3_dott, V_p_xt, V_p_yt, V_p_zt, w_p_xt, w_p_yt, w_p_zt];
 
ls =     [l1,         l2,     l3];
l_dots = [l1_dot, l2_dot, l3_dot];

thetas     = [    theta1,     theta2,     theta3];
theta_dots = [theta1_dot, theta2_dot, theta3_dot];

Q = [T1 T2 F3 F4 F5 0 0 0];

% w = [wx wy wz].';

%% Platform Only (No other joints) Forward Computation Calcs
% This gets the parameters for wx_platform, wy_platform, wz_platform 

fprintf("Calculating forward kinematics\n");

for i = 1:3
    B_temp{i} = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(ls(i));
    B{i} = B_temp{i}(1:3,4);
    b_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    b{i} = b_temp{i}(1:3,4);
end

% P_c = (B{1} + B{2} + B{3})/3;

% for i = 1:2
%     eqs_(-2+3*i:0+3*i) = T(1:3,1:3)*b{i} + P_c == B{i};
% end
% 
% [A_lin,B_lin] = equationsToMatrix(eqs_, [a1 a2 a3 o1 o2 o3]);
% 
% X = linsolve(A_lin,B_lin);
% 
% X(1:3) = X(1:3)/norm(X(1:3));
% X(4:6) = X(4:6)/norm(X(4:6));
% 
% a1_ = simplify(X(1));
% a2_ = simplify(X(2));
% a3_ = simplify(X(3));
% o1_ = simplify(X(4));
% o2_ = simplify(X(5));
% o3_ = simplify(X(6));
% n1_ = simplify( o2_*a3_ - a2_*o3_);
% n2_ = simplify(-o1_*a3_ + a1_*o3_);
% n3_ = simplify( o1_*a2_ - o2_*a1_);
% 
% unit_vecs = [a1_ a2_ a3_ o1_ o2_ o3_ n1_ n2_ n3_];

%%

fprintf("Calculating w_platform\n");

% px = P_c(1); py = P_c(2); pz = P_c(3);
% 
% V_c = [0 0 0].';
% 
% for i = 1:length(qs)
%     V_c = V_c + diff(P_c,qs(i))*q_dots(i);
% end
% 
% V_c = simplify(V_c);
% 
% for i = 1:3
%     r_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
%     r_{i} = r_temp{i}(1:3,4);
%     R_0_2{i} = Rx(alpha_5 + gamma(i))*Rz(-pi/2+thetas(i));
%     phi_l{i} = R_0_2{i}(1:3,1);
%     phi_z{i} = R_0_2{i}(1:3,3);
% end
% 
% for i = 1:3
%    Vbi_1{i} = Tinv(1:3,1:3)*V_c + cross(w,r_{i});
%    Vbi_2{i} = Tinv(1:3,1:3)*(l_dots(i)*phi_l{i} + cross(theta_dots(i)*phi_z{i},ls(i)*phi_l{i}));
%    eq{i} = Vbi_1{i} - Vbi_2{i} == 0;
% end
% 
% wx_solved = simplify(solve(eq{1}(2),wx));
% [Alin,Blin] = equationsToMatrix([eq{1}(1);eq{2}(1)], [wy wz]);
% wy_wz = Alin\Blin;
% wy_solved = simplify(wy_wz(1));
% wz_solved = simplify(wy_wz(2));
% 
% % substitue in unit vectors to get into a form we can use for dynamics
% wx_solved = simplify(subs(wx_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
% wy_solved = simplify(subs(wy_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
% wz_solved = simplify(subs(wz_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));

%% Position and Velocity of the Elbow including other links

fprintf("Calculating p_elbow, v_elbow given other joints\n");

E_com_in_frame = [-0.07272380,-0.00932668,0.18779160];

P_e_com_temp = Rz(qe)*TRANSx(E_com_in_frame(1))*TRANSy(E_com_in_frame(2))*TRANSz(E_com_in_frame(3));
P_e_com = simplify(P_e_com_temp(1:3,4));

V_e_com = [0 0 0].';

for i = 1:length(qs)
    V_e_com = V_e_com + diff(P_e_com,qs(i))*q_dots(i);
end

V_e_com = simplify(V_e_com);

w_e = [0,0,qe_dot].';

%% Position and Velocity of the Forearm including other links

fprintf("Calculating p_forearm, v_forearm given other joints\n");

F_com_in_frame = [-0.00401306,0,0]; % Could add y, and z. Keeping 0 for simplification for now bc they are small

P_f_com_temp = Rz(qe)*TRANSx(a4)*Rx(qf)*TRANSx(F_com_in_frame(1));
P_f_com = simplify(P_f_com_temp(1:3,4));

V_f_com = [0 0 0].';

for i = 1:length(qs)
    V_f_com = V_f_com + diff(P_f_com,qs(i))*q_dots(i);
end

V_f_com = simplify(V_f_com);

R_e_f_temp = Rx(qf);

w_f = [qf_dot,0,0].' + R_e_f_temp(1:3,1:3)*w_e;

%% Position and Velocity of the Linear Links including other links

fprintf("Calculating p_links, v_links given other joints\n");

for i = 1:3
    d{i} = (ls(i) - l_offset);
    P_l_com_temp{i} = Rz(qe)*TRANSx(a4)*Rx(qf + alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(d{i});
    P_l_com{i} = P_l_com_temp{i}(1:3,4);
end

V_l_com{1} = [0 0 0].';
V_l_com{2} = [0 0 0].';
V_l_com{3} = [0 0 0].';

for i = 1:3
    for j = 1:length(qs)
        V_l_com{i} = V_l_com{i} + diff(P_l_com{i},qs(j))*q_dots(j);
    end
    V_l_com{i} = simplify(V_l_com{i});
end

for i = 1:3
    R_f_l_temp{i} = Rx(alpha_5+gamma(i))*Rz(-pi/2+qs(i+5));
    w_l{i} = [0,0,q_dots(i+5)].' + R_f_l_temp{i}(1:3,1:3)*w_f;
end

%% Position and Velocity of the Platform including other links

fprintf("Calculating p_platform, v_platform given other joints\n");

% for i = 1:3
%     Bcom_temp{i} = Rz(qe)*TRANSx(a4)*Rx(qf + alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(ls(i));
%     Bcom{i} = Bcom_temp{i}(1:3,4);
% end
% 
% P_p_com = simplify((Bcom{1} + Bcom{2} + Bcom{3})/3);
% 
% V_p_com = [0 0 0].';
% 
% for i = 1:length(qs)
%     V_p_com = V_p_com + diff(P_p_com,qs(i))*q_dots(i);
% end
% 
% V_p_com = simplify(V_p_com);
% 
% R_f_p = [n1_, n2_, n3_; 
%          o1_, o2_, o3_; 
%          a1_, a2_, a3_];

% w_p = R_f_p*w_f + [wx_solved,wy_solved,wz_solved].';

P_f_o = Rz(qe)*TRANSx(a4)*Rx(qf);
P_f_o = simplify(P_f_o(1:3,4));

V_f_o = [0 0 0].';

%velocity of forearm origin instead of the COM
for i = 1:length(qs)
    V_f_o = V_f_o + diff(P_f_o,qs(i))*q_dots(i);
end

% rotation from origin to the forearm
R_origin_f = Rz(qe)*Rx(qf);
P_p_com_temp = Rz(qe)*TRANSx(a4)*Rx(qf);

P_p_com = P_p_com_temp(1:3,1:3)*[P_p_x, P_p_y, P_p_z].';
P_p_com = P_p_com(1:3);
V_p_com = V_f_o + R_origin_f(1:3,1:3)*[V_p_x, V_p_y, V_p_z].';

R_p_com = [R_p_x, R_p_y, R_p_z];

R_f_p = Ry(R_p_y)*Rz(R_p_z)*Rx(R_p_x);

w_p = R_f_p(1:3,1:3)*w_f + [w_p_x, w_p_y, w_p_z].';

%%
fprintf("Calculating Lagrangian\n");

T = 1/2*Mp*(V_p_com(1)^2 + V_p_com(2)^2 + V_p_com(3)^2)...
  + 1/2*(Ixx_p*w_p(1)^2 + Iyy_p*w_p(2)^2 + Izz_p*w_p(3)^2)...
  + 1/2*Mf*(V_f_com(1)^2 + V_f_com(2)^2 + V_f_com(3)^2)...
  + 1/2*(Ixx_f*w_f(1)^2 + Iyy_f*w_f(2)^2 + Izz_f*w_f(3)^2)...
  + 1/2*Me*(V_e_com(1)^2 + V_e_com(2)^2 + V_e_com(3)^2)...
  + 1/2*(Ixx_e*w_e(1)^2 + Iyy_e*w_e(2)^2 + Izz_e*w_e(3)^2)...
  + 1/2*Ml*(V_l_com{1}(1)^2 + V_l_com{1}(2)^2 + V_l_com{1}(3)^2)...
  + 1/2*(Ixx_l*w_l{1}(1)^2 + Iyy_l*w_l{1}(2)^2 + Izz_l*w_l{1}(3)^2)...
  + 1/2*Ml*(V_l_com{2}(1)^2 + V_l_com{2}(2)^2 + V_l_com{2}(3)^2)...
  + 1/2*(Ixx_l*w_l{2}(1)^2 + Iyy_l*w_l{2}(2)^2 + Izz_l*w_l{2}(3)^2)...
  + 1/2*Ml*(V_l_com{3}(1)^2 + V_l_com{3}(2)^2 + V_l_com{3}(3)^2)...
  + 1/2*(Ixx_l*w_l{3}(1)^2 + Iyy_l*w_l{3}(2)^2 + Izz_l*w_l{3}(3)^2);

P = g * (Mp*P_p_com(2) + Mf*P_f_com(2) + Me*P_e_com(2) + Ml*(P_l_com{1}(2)+P_l_com{2}(2)+P_l_com{3}(2)));

L = T-P;
% L = subs(L,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs);

%% COMPUTING EOMS

fprintf("Computing EOMs\n");

for i = 1:length(qs)
    dldqdot{i} = diff(L,q_dots(i));
    dldqdot_t{i} = subs(dldqdot{i},[qs q_dots], [qst q_dotst]);
    d_dldqdot_t{i} = diff(dldqdot_t{i},t);
    d_dldqdot{i} = subs(d_dldqdot_t{i},[qst q_dotst diff([qst q_dotst],t)],[qs q_dots q_dots q_d_dots]);
    
    EOMS(i) = d_dldqdot{i} - diff(L,qs(i));
end

[M, V, G] = separate_mvg_no_simp(EOMS.',q_d_dots.',g);
%%

fprintf("Converting V vector to a square matrix\n");

Vsquare = sym(zeros(length(qs),length(qs)));

for i = 1:length(qs)
    test = V(i);
    for j = 1:length(qs)
        [coe,terms] = coeffs(test,q_dots(j));
        for k = 1:length(terms)
            if subs(terms(k),q_dots(j),0) == 0
                Vsquare(i,j) = Vsquare(i,j) + coe(k)*simplify(terms(k)/q_dots(j));
            elseif subs(terms(k),q_dots(j),0) ~= 0
                test = coe(k);
            end
        end
    end
end

%%
fprintf("Computing and incorporating Rho\n");

% fk(1,1) = simplify(expand(sqrt((B{1}(1)-B{2}(1))^2+(B{1}(2)-B{2}(2))^2+(B{1}(3)-B{2}(3))^2)-sqrt(3)*r));
% fk(2,1) = simplify(expand(sqrt((B{3}(1)-B{2}(1))^2+(B{3}(2)-B{2}(2))^2+(B{3}(3)-B{2}(3))^2)-sqrt(3)*r));
% fk(3,1) = simplify(expand(sqrt((B{3}(1)-B{1}(1))^2+(B{3}(2)-B{1}(2))^2+(B{3}(3)-B{1}(3))^2)-sqrt(3)*r));

for i = 1:3
    B_temp{i} = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(ls(i));
    B{i} = B_temp{i}(1:3,4);
    b_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    b{i} = b_temp{i}(1:3,4);
    
    fk(3*(i-1)+1:3*(i-1)+3,1) = B{i} - R_f_p(1:3,1:3)*b{i} - [P_p_x, P_p_y, P_p_z].';
end

psi = [fk;qs(1:(length(qs)-9)).'];
psi_dq = simplify(jacobian(psi,qs));
psi_dq_t = subs(psi_dq,qs,qst);
psi_dq_dt = subs(diff(psi_dq_t,t),[qst diff(qst,t)],[qs, q_dots]);


% rho = inv(psi_dq)*[zeros(9,(length(qs)-9));eye((length(qs)-9),(length(qs)-9))];
% rho_t = subs(rho,qs,qst);
% rho_dt = subs(diff(rho_t,t),[qst diff(qst,t)],[qs, q_dots]);

% Mpar = rho.'*M*rho;
% Vpar = rho.'*Vsquare*rho + rho.'*M*rho_dt;
% Gpar = rho.'*G;

%% Substitute in actual values
fprintf("Substituting in values\n")

% Variables
l_offset = 0.08803884;
r        = 0.05288174521;
R        = 0.1044956;
alpha_5  = 0.094516665; 
alpha_13 = 5*pi/180;
a56      = 3.8340e-04;
a4       = 0.159385;
g        = 9.81;

% Elbow
Ixx_e = 0.09961102;
Iyy_e = 0.13932454;
Izz_e = 0.06204712;
Me    = 7.10726178;

% Forearm
Ixx_f = 0.01746898;
Iyy_f = 0.00902027;
Izz_f = 0.00893149;
Mf    = 1.45613142;

% Prismatic Links
Ixx_l = 0.00000664;
Iyy_l = 0.00029962;
Izz_l = 0.00029713;
Ml    = 0.14820004;

% Platform
Ixx_p = 0.00090728;
Iyy_p = 0.00048834;
Izz_p = 0.00107519;
Mp    = 0.36486131;

% Substitute in real values
% Mpar    = subs(Mpar);
% Vpar    = subs(Vpar);
% Gpar    = subs(Gpar);
M         = subs(M);
Vsquare   = subs(Vsquare);
G         = subs(G);
psi       = subs(psi);
psi_dq    = subs(psi_dq);
psi_dq_dt = subs(psi_dq_dt);

% rho     = subs(rho);
% rho_dt  = subs(rho_dt);

%%
% mass_props = [Ixxf Ixx Iyy Izz Mp Ml l_offset g];
% mass_values = [0.01747099 0.00091827 0.00080377 0.00138661 0.36486131 0.14820004 0.08803884 9.81];
q_vals =      [0.0, 0.0, 0.1, 0.1, 0.1, 1.02844, 1.02844, 1.02844, .0856, 0,    0, 0,     0, 0];
q_dot_vals  = [0.2, 0.1,   0,   0,   0,   -0.01,       0,       0,     0, 0, 0.02, 0, -0.01, 0];
q_ddot_vals = [  0,   0,   0,   0,   0,       0,       0,       0,     0, 0,    0, 0,     0, 0];
test = subs(M,[qs q_dots q_d_dots], [q_vals q_dot_vals q_ddot_vals]);
% psi_dq_subs = subs(psi_dq,[qs q_dots q_d_dots], [q_vals q_dot_vals q_ddot_vals]);
% psi_dq_dt_subs = subs(psi_dq_dt,[qs q_dots q_d_dots], [q_vals q_dot_vals q_ddot_vals]);
% test2 = -inv(psi_dq_subs)*psi_dq_dt_subs*inv(psi_dq_subs)*[zeros(9,(length(qs)-9));eye((length(qs)-9),(length(qs)-9))];

vpa(test,4)
% vpa(test2,4)

%% Write parallel structure to file

fprintf("Writing to file\n");

M_write = ccode(vpa(M,5));
fileID = fopen("DynamicEqs/M.txt","w");
fprintf(fileID,M_write);

V_write = ccode(vpa(Vsquare,5));
fileID = fopen("DynamicEqs/V.txt","w");
fprintf(fileID,V_write);

G_write = ccode(vpa(G,5));
fileID = fopen("DynamicEqs/G.txt","w");
fprintf(fileID,G_write);

psi_write = ccode(vpa(psi,5));
fileID = fopen("DynamicEqs/psi.txt","w");
fprintf(fileID,psi_write);

psi_dq_write = ccode(vpa(psi_dq,5));
fileID = fopen("DynamicEqs/psi_dq.txt","w");
fprintf(fileID,psi_dq_write );

psi_dq_dt_write = ccode(vpa(psi_dq_dt,5));
fileID = fopen("DynamicEqs/psi_dq_dt.txt","w");
fprintf(fileID,psi_dq_dt_write);