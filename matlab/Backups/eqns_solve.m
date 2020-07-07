function F = eqns_solve(theta, l)

R = 0.1044956; % [m]
r = 0.05288174521; % [m]
a5 = 0.0268986; % [m]
a6 = 0.0272820; % [m]
a56 = -(a5-a6); % [m]

t1 = theta(1);
t2 = theta(2);
t3 = theta(3);

l1 = l(1);
l2 = l(2);
l3 = l(3);

F(1) = (3*R^2 + 3*a56^2 + l1^2 + l2^2 - 3*R*l1*sin(t1) - 3*R*l2*sin(t2) - (l1*l2*cos(t1 - t2))/2 - (3*l1*l2*cos(t1 + t2))/2 - 3^(1/2)*a56*l1*sin(t1) + 3^(1/2)*a56*l2*sin(t2))^(1/2) - 3^(1/2)*r;
F(2) = (3*R^2 + 3*a56^2 + l2^2 + l3^2 - 3*R*l2*sin(t2) - 3*R*l3*sin(t3) - (l2*l3*cos(t2 - t3))/2 - (3*l2*l3*cos(t2 + t3))/2 - 3^(1/2)*a56*l2*sin(t2) + 3^(1/2)*a56*l3*sin(t3))^(1/2) - 3^(1/2)*r;
F(3) = (3*R^2 + 3*a56^2 + l1^2 + l3^2 - 3*R*l1*sin(t1) - 3*R*l3*sin(t3) - (l1*l3*cos(t1 - t3))/2 - (3*l1*l3*cos(t1 + t3))/2 + 3^(1/2)*a56*l1*sin(t1) - 3^(1/2)*a56*l3*sin(t3))^(1/2) - 3^(1/2)*r;
end

    
    