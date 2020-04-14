function [T] = Rx(theta)

% returns Homogeneous Transform about x given theta

T = [1 0           0           0;
      0 cos(theta) -sin(theta) 0;
      0 sin(theta) cos(theta)  0;
      0 0           0          1];
  
    T = zeromat(T,1e-10);
end