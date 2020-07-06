function [T] = Rz(theta)

% returns Homogeneous Transform about z given theta

    T = [ cos(theta)  -sin(theta) 0  0;
          sin(theta)  cos(theta)  0  0;
          0           0           1  0;
          0           0           0  1];

    T = zeromat(T,1e-10);
end