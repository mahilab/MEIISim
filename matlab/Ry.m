function [T] = Ry(theta)

% returns Homogeneous Transform about y given theta

    T = [ cos(theta)  0 sin(theta) 0;
          0           1 0          0;
          -sin(theta) 0 cos(theta) 0;
          0           0 0          1];

  T = zeromat(T,1e-10);
end