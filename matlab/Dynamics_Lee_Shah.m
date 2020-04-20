syms l1_dot l2_dot l3_dot theta1_dot theta2_dot theta3_dot

P_c = [(l1*cos(theta1))/3 + (l2*cos(theta2))/3 + (l3*cos(theta3))/3;(R*cos(alpha_5))/3 + (a56*sin(alpha_5))/3 - (R*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 - (R*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (a56*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 - (l1*cos(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (l3*sin(theta3)*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3;(R*sin(alpha_5))/3 - (a56*cos(alpha_5))/3 - (R*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3 - (R*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 - (3^(1/2)*sin(alpha_5))/2))/3 + (a56*(cos(alpha_5)/2 + (3^(1/2)*sin(alpha_5))/2))/3 - (l1*sin(alpha_5)*sin(theta1))/3 + (l2*sin(theta2)*(sin(alpha_5)/2 + (3^(1/2)*cos(alpha_5))/2))/3 + (l3*sin(theta3)*(sin(alpha_5)/2 - (3^(1/2)*cos(alpha_5))/2))/3];

V_c = [0 0 0].'



