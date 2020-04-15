function [ new_point ] = rotate_point_vector( point, vector, center, angle)
%ROTATE_POINT_VECTOR Summary of this function goes here
%   Detailed explanation goes here
    a = vector(1); b = vector(2); c = vector(3);
    x1 = center(1); y1 = center(2); z1 = center(3);
    d = sqrt(b^2+c^2);
    T = [1 0 0 -x1;
         0 1 0 -y1;
         0 0 1 -z1;
         0 0 0  1];
    Rx = [1 0 0 0;
          0 c/d -b/d 0;
          0 b/d c/d 0;
          0 0 0 1];
    Ry = [d 0 -a 0;
          0 1 0 0;
          a 0 d 0;
          0 0 0 1];
    Rz = [cos(angle), -sin(angle), 0, 0;
          sin(angle), cos(angle), 0, 0;
          0 0 1 0;
          0 0 0 1];
    new_point = inv(T)*inv(Rx)*inv(Ry)*Rz*Rz*Rx*T*[point.';1];
    new_point = new_point(1:3);
end

