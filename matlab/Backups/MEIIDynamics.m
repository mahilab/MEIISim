clc; clear all; close all;

syms tau1 eta1 Jm1 q1 q1d q1dd mass1 b1 fk1 fenv1x fenv1y fenv1z menv1x menv1y menv1z... % Link 1
     tau2 eta2 Jm2 q2 q2d q2dd mass2 b2 fk2 fenv2x fenv2y fenv2z menv2x menv2y menv2z... % Link 2
     tau3 eta3 Jm3 q3 q3d q3dd mass3 b3 fk3 fenv3x fenv3y fenv3z menv3x menv3y menv3z... % Link 3
     tau4 eta4 Jm4 q4 q4d q4dd mass4 b1 fk1 fenv4x fenv4y fenv4z menv4x menv4y menv4z... % Link 4
     tau5 eta5 Jm5 q5 q5d q5dd mass5 b2 fk2 fenv5x fenv5y fenv5z menv5x menv5y menv5z... % Link 5
     tau6 eta6 Jm6 q6 q6d q6dd mass6 b3 fk3 fenv6x fenv6y fenv6z menv6x menv6y menv6z... % Link 6
     tau7 eta7 Jm7 q7 q7d q7dd mass7 b1 fk1 fenv7x fenv7y fenv7z menv7x menv7y menv7z... % Link 7
     tau8 eta8 Jm8 q8 q8d q8dd mass8 b1 fk1 fenv8x fenv8y fenv8z menv8x menv8y menv8z... % Link 8
     r1x r1y r1z ms1x ms1y ms1z Ic1xx Ic1xy Ic1xz Ic1yy Ic1yz Ic1zz ... % Link 1
     r2x r2y r2z ms2x ms2y ms2z Ic2xx Ic2xy Ic2xz Ic2yy Ic2yz Ic2zz ... % Link 2
     r3x r3y r3z ms3x ms3y ms3z Ic3xx Ic3xy Ic3xz Ic3yy Ic3yz Ic3zz ... % Link 3
     r4x r4y r2z ms4x ms4y ms4z Ic4xx Ic4xy Ic4xz Ic4yy Ic4yz Ic4zz ... % Link 4
     r5x r5y r3z ms5x ms5y ms5z Ic5xx Ic5xy Ic5xz Ic5yy Ic5yz Ic5zz ... % Link 5
     r6x r6y r2z ms6x ms6y ms6z Ic6xx Ic6xy Ic6xz Ic6yy Ic6yz Ic6zz ... % Link 6
     r7x r7y r3z ms7x ms7y ms7z Ic7xx Ic7xy Ic7xz Ic7yy Ic7yz Ic7zz ... % Link 7
     r8x r8y r2z ms8x ms8y ms8z Ic8xx Ic8xy Ic8xz Ic8yy Ic8yz Ic8zz ... % Link 8
     R a4 a56 alpha_5... % Parameters for modified dh parameters
     g

% Initializing matrices of symbolic parameters 
ms1 = [ms1x ms1y ms1z].';
ms2 = [ms2x ms2y ms2z].';
ms3 = [ms3x ms3y ms3z].';
ms4 = [ms4x ms4y ms4z].';
ms5 = [ms5x ms5y ms5z].';
ms6 = [ms6x ms6y ms6z].';
ms7 = [ms7x ms7y ms7z].';
ms8 = [ms8x ms8y ms8z].';

fenv1 = [fenv1x fenv1y fenv1z].';
fenv2 = [fenv2x fenv2y fenv2z].';
fenv3 = [fenv3x fenv3y fenv3z].';
fenv4 = [fenv4x fenv4y fenv4z].';
fenv5 = [fenv5x fenv5y fenv5z].';
fenv6 = [fenv6x fenv6y fenv6z].';
fenv7 = [fenv7x fenv7y fenv7z].';
fenv8 = [fenv8x fenv8y fenv8z].';

menv1 = [menv1x menv1y menv1z].';
menv2 = [menv2x menv2y menv2z].';
menv3 = [menv3x menv3y menv3z].';
menv4 = [menv4x menv4y menv4z].';
menv5 = [menv5x menv5y menv5z].';
menv6 = [menv6x menv6y menv6z].';
menv7 = [menv7x menv7y menv7z].';
menv8 = [menv8x menv8y menv8z].';

Tau  = [tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8].'; % Vector of torques provided at each coordinate
Eta  = [eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8].'; % Vector of transmission ratios at each coordinate
Jm   = [Jm1 Jm2 Jm3 Jm4 Jm5 Jm6 Jm7 Jm8].';         % Vector of rotor inertias for each joint
mass = [mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8].'; % Mass of each link

q    = [q1 q2 q3 q4 q5 q6 q7 q8].';                 % movement variable at each joint
q_d  = [q1d q2d q3d q4d q5d q6d q7d q8d].';         % first derivative of movement variable
q_dd = [q1dd q2dd q3dd q4dd q5dd q6dd q7dd q8dd].'; % second derivative of movement variable

B = [b1 b2 b3].';
Fk = [fk1 fk2 fk3].';

Ic1 = [Ic1xx -Ic1xy -Ic1xz;
    -Ic1xy Ic1yy -Ic1yz;
    -Ic1xz -Ic1yz Ic1xx];

Ic2 = [Ic2xx -Ic2xy -Ic2xz;
    -Ic2xy Ic2yy -Ic2yz;
    -Ic2xz -Ic2yz Ic2xx];

Ic3 = [Ic3xx -Ic3xy -Ic3xz;
    -Ic3xy Ic3yy -Ic3yz;
    -Ic3xz -Ic3yz Ic3xx];
    
Ic4 = [Ic4xx -Ic4xy -Ic4xz;
    -Ic4xy Ic4yy -Ic4yz;
    -Ic4xz -Ic4yz Ic4xx];

Ic5 = [Ic5xx -Ic5xy -Ic5xz;
    -Ic5xy Ic5yy -Ic5yz;
    -Ic5xz -Ic5yz Ic5xx];

Ic6 = [Ic6xx -Ic6xy -Ic6xz;
    -Ic6xy Ic6yy -Ic6yz;
    -Ic6xz -Ic6yz Ic6xx];

Ic7 = [Ic7xx -Ic7xy -Ic7xz;
    -Ic7xy Ic7yy -Ic7yz;
    -Ic7xz -Ic7yz Ic7xx];

Ic8 = [Ic8xx -Ic8xy -Ic8xz;
    -Ic8xy Ic8yy -Ic8yz;
    -Ic8xz -Ic8yz Ic8xx];

ms = {ms1 ms2 ms3 ms4 ms5 ms6 ms7 ms8};
Ic = {Ic1 Ic2 Ic3 Ic4 Ic5 Ic6 Ic7 Ic8};
menv = {menv1 menv2 menv3 menv4 menv5 menv6 menv7 menv8};
fenv = {fenv1 fenv2 fenv3 fenv4 fenv5 fenv6 fenv7 fenv8};

%% Forward kinematics

% Modified DH Paramaters for each link
% Link   1   2           3      4     5       6     7       8
a     = [0,  1,          2,     3,    2,      5,    2,      7   ];
mu    = [1,  1,          0,     1,    0,      1,    0,      1   ];
sigma = [0,  0,          0,     1,    0,      1,    0,      1   ];
gamma = [0,  0,          0,     0,    2*pi/3, 0,    4*pi/3, 0   ];
b     = [0,  0,          0,     0,    0,      0,    0,      0   ];
alpha = [0,  pi/2,       -pi/2, pi/2, -pi/2,  pi/2, -pi/2,  pi/2];
d     = [0,  0,          R,     0,    R,      0,    R,      0   ];
theta = [q1, q2+alpha_5, q3,    0,    q5,     0,    q7,     0   ];
r     = [0,  a4,         a56,   q4,   a56,    q6,   a56,    q8  ];
 
num_links = 8;
 
T = cell(num_links,1);

for i=1:num_links
    T{i} = mdhtf(gamma(i), b(i), alpha(i), d(i), theta(i), r(i));
end

% g = 9.80665;
g0 = [-g 0 0].';

%% Run dynamics
MVG = newton_euler_dynamics(a, sigma, T, mass, ms, Ic, fenv, menv, q_d, q_dd, g0);
MVG = MVG(2:9);
MVG = simplify(expand(MVG));

%% Separate MVG into M, V, and G
[M,V,G] = separate_mvg(MVG,q_dd,g);


