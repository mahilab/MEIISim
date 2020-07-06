clear all;

syms qf l1 l2 l3 theta1 theta2 theta3...
     qf_dot l1_dot l2_dot l3_dot theta1_dot theta2_dot theta3_dot...
     qf_d_dot l1_d_dot l2_d_dot l3_d_dot theta1_d_dot theta2_d_dot theta3_d_dot...
     qft(t) l1t(t) l2t(t) l3t(t) theta1t(t) theta2t(t) theta3t(t)...
     qf_dott(t) l1_dott(t) l2_dott(t) l3_dott(t) theta1_dott(t) theta2_dott(t) theta3_dott(t)...
     wx wy wz vx vy vz px py pz...
     r R alpha_5 a56 alpha_13...
     T1 F1 F2 F3...
     n1 n2 n3 o1 o2 o3 a1 a2 a3...
     g
 
T = [n1 o1 a1 px;
     n2 o2 a2 py;
     n3 o3 a3 pz;
      0  0  0  1];
 
Tinv = [T(1:3,1:3).',-T(1:3,1:3).'*T(1:3,4);T(4,:)];
 
gamma = [0, -2*pi/3, 2*pi/3];

qs       = [      qf,       l1,       l2,       l3,       theta1,       theta2,       theta3];
q_dots   = [  qf_dot,   l1_dot,   l2_dot,   l3_dot,   theta1_dot,   theta2_dot,   theta3_dot];
q_d_dots = [qf_d_dot, l1_d_dot, l2_d_dot, l3_d_dot, theta1_d_dot, theta2_d_dot, theta3_d_dot];

qst     = [    qft,     l1t,     l2t,     l3t,     theta1t,     theta2t,     theta3t];
q_dotst = [qf_dott, l1_dott, l2_dott, l3_dott, theta1_dott, theta2_dott, theta3_dott];
 
ls =     [l1,         l2,     l3];
l_dots = [l1_dot, l2_dot, l3_dot];

thetas     = [    theta1,     theta2,     theta3];
theta_dots = [theta1_dot, theta2_dot, theta3_dot];

Q = [T1 F1 F2 F3 0 0 0];

w = [wx wy wz].';

%% Platform Only (No other joints) Forward Computation Calcs
% This gets the parameters for wx_platform, wy_platform, wz_platform 

fprintf("Calculating forward kinematics\n");

for i = 1:3
    B_temp{i} = Rx(alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(ls(i));
    B{i} = B_temp{i}(1:3,4);
    b_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    b{i} = b_temp{i}(1:3,4);
end

P_c = (B{1} + B{2} + B{3})/3;

for i = 1:2
    eqs_(-2+3*i:0+3*i) = T(1:3,1:3)*b{i} + P_c == B{i};
end

[A_lin,B_lin] = equationsToMatrix(eqs_, [a1 a2 a3 o1 o2 o3]);

X = linsolve(A_lin,B_lin);

X(1:3) = X(1:3)/norm(X(1:3));
X(4:6) = X(4:6)/norm(X(4:6));

a1_ = simplify(X(1));
a2_ = simplify(X(2));
a3_ = simplify(X(3));
o1_ = simplify(X(4));
o2_ = simplify(X(5));
o3_ = simplify(X(6));
n1_ = simplify( o2_*a3_ - a2_*o3_);
n2_ = simplify(-o1_*a3_ + a1_*o3_);
n3_ = simplify( o1_*a2_ - o2_*a1_);

unit_vecs = [a1_ a2_ a3_ o1_ o2_ o3_ n1_ n2_ n3_];

%%

fprintf("Calculating w_platform\n");

px = P_c(1); py = P_c(2); pz = P_c(3);

V_c = [0 0 0].';

for i = 1:length(qs)
    V_c = V_c + diff(P_c,qs(i))*q_dots(i);
end

V_c = simplify(V_c);

for i = 1:3
    r_temp{i} = Rx(alpha_13 + gamma(i))*TRANSy(r);
    r_{i} = r_temp{i}(1:3,4);
    R_0_2{i} = Rx(alpha_5 + gamma(i))*Rz(-pi/2+thetas(i));
    phi_l{i} = R_0_2{i}(1:3,1);
    phi_z{i} = R_0_2{i}(1:3,3);
end

for i = 1:3
   Vbi_1{i} = Tinv(1:3,1:3)*V_c + cross(w,r_{i});
   Vbi_2{i} = Tinv(1:3,1:3)*(l_dots(i)*phi_l{i} + cross(theta_dots(i)*phi_z{i},ls(i)*phi_l{i}));
   eq{i} = Vbi_1{i} - Vbi_2{i} == 0;
end

wx_solved = simplify(solve(eq{1}(2),wx));
[Alin,Blin] = equationsToMatrix([eq{1}(1);eq{2}(1)], [wy wz]);
wy_wz = Alin\Blin;
wy_solved = simplify(wy_wz(1));
wz_solved = simplify(wy_wz(2));

% substitue in unit vectors to get into a form we can use for dynamics
wx_solved = simplify(subs(wx_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
wy_solved = simplify(subs(wy_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));
wz_solved = simplify(subs(wz_solved,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs));

%% Position and Velocity of the Platform including other links

fprintf("Calculating p_platform, v_platform given other joints\n");

for i = 1:3
    Bcom_temp{i} = Rx(qf + alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(ls(i));
    Bcom{i} = Bcom_temp{i}(1:3,4);
end

P_com = simplify((Bcom{1} + Bcom{2} + Bcom{3})/3);

V_com = [0 0 0].';

for i = 1:length(qs)
    V_com = V_com + diff(P_com,qs(i))*q_dots(i);
end

V_com = simplify(V_com);

%%
fprintf("Calculating Lagrangian\n");

syms Ixx Iyy Izz Mp Ml l_offset Ixxf Iyyf Izzf

sum_By  = 0;

V_l{1} = [0;0;0];
V_l{2} = [0;0;0];
V_l{3} = [0;0;0];

for i = 1:3
    d{i} = (ls(i) - l_offset);
    B_temp_y{i} = Rx(qf + alpha_5 + gamma(i))*TRANSy(R)*TRANSz(-a56)*Rz(-pi/2+thetas(i))*TRANSx(d{i});
    By{i} = B_temp_y{i}(2,4);
    sum_By = sum_By + By{i};
    
    for j= 1:length(qs)
        V_l{i} = V_l{i} + diff(B_temp_y{i}(1:3,4),qs(j))*q_dots(j);
    end
end

T = 1/2*Mp*(V_com(1)^2 + V_com(2)^2 + V_com(3)^2)...
  + 1/2*(Ixx*(wx_solved+qf_dot)^2 + Iyy*wy_solved^2 + Izz*wz_solved^2)...
  + Ixxf*(qf_dot^2)...
  + 1/2*Ml*(V_l{1}(1)^2 + V_l{1}(2)^2 + V_l{1}(3)^2)...
  + 1/2*Ml*(V_l{2}(1)^2 + V_l{2}(2)^2 + V_l{2}(3)^2)...
  + 1/2*Ml*(V_l{3}(1)^2 + V_l{3}(2)^2 + V_l{3}(3)^2);

P = Mp*g*py + simplify(Ml*g*sum_By);

L = T-P;
L = subs(L,[a1 a2 a3 o1 o2 o3 n1 n2 n3],unit_vecs);

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

fk(1,1) = simplify(expand(sqrt((B{1}(1)-B{2}(1))^2+(B{1}(2)-B{2}(2))^2+(B{1}(3)-B{2}(3))^2)-sqrt(3)*r));
fk(2,1) = simplify(expand(sqrt((B{3}(1)-B{2}(1))^2+(B{3}(2)-B{2}(2))^2+(B{3}(3)-B{2}(3))^2)-sqrt(3)*r));
fk(3,1) = simplify(expand(sqrt((B{3}(1)-B{1}(1))^2+(B{3}(2)-B{1}(2))^2+(B{3}(3)-B{1}(3))^2)-sqrt(3)*r));

psi = [fk;qs(1:(length(qs)-3)).'];
psi_dq = simplify(jacobian(psi,qs));
rho = simplify(inv(psi_dq)*[zeros(3,(length(qs)-3));eye((length(qs)-3),(length(qs)-3))]);
rho_t = subs(rho,qs,qst);
rho_dt = subs(diff(rho_t,t),[qst diff(qst,t)],[qs, q_dots]);

Mpar = rho.'*M*rho;
Vpar = rho.'*Vsquare*rho + rho.'*M*rho_dt;
Gpar = rho.'*G;

%% Substitute in actual values
fprintf("Substituting in values\n")

Ixx = 0.00091827;
Iyy = 0.00080377;
Izz = 0.00138661;
Mp = 0.36486131;
Ml = 0.14820004;
l_offset = 0.08803884;
r        = 0.05288174521;
R        = 0.1044956;
alpha_5  = 0.094516665; 
alpha_13 = 5*pi/180;
a56      = 3.8340e-04;
g        = 9.81;

Ixxf = 0.01747099;
Iyyf = 0.00904374;
Izzf = 0.00895693;

Mpar = subs(Mpar);
Vpar = subs(Vpar);
Gpar = subs(Gpar);

%%
% mass_props = [Ixxf Ixx Iyy Izz Mp Ml l_offset g];
% mass_values = [0.01747099 0.00091827 0.00080377 0.00138661 0.36486131 0.14820004 0.08803884 9.81];
% q_q_dot_vals = [0.0, 0.1, 0.1, 0.1, 1.0284, 1.0284, 1.0284, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
% q_ddot_vals = [0.00 0.00 0.00 0.00 0.00 0.00 0.00];
% variables = [r R alpha_5 alpha_13 a56];
% values = [0.05288174521 0.1044956 0.094516665 5*pi/180 3.8340e-04];
% test = subs(subs(Gpar),[qs q_dots q_d_dots variables mass_props], [q_q_dot_vals q_ddot_vals values mass_values]);
% 
% vpa(test,4)

%% Write parallel structure to file

fprintf("Writing to file\n");

for i = 1:(length(qs)-3)
    for j = 1:(length(qs)-3)
        M_write = ccode(vpa(Mpar(i,j),3));
        fileID = fopen("DynamicEqs/Mpar" + int2str(i) + int2str(j) + ".txt","w");
        fprintf(fileID,M_write);
        V_write = ccode(vpa(Vpar(i,j),3));
        fileID = fopen("DynamicEqs/Vpar" + int2str(i) + int2str(j) + ".txt","w");
        fprintf(fileID,V_write);
    end
    G_write = ccode(vpa(Gpar(i),3));
    fileID = fopen("DynamicEqs/Gpar" + int2str(i) + ".txt","w");
    fprintf(fileID,G_write);
end