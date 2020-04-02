function [ T ] = mdhtf( gamma, b, alpha, d, theta, r)
%mdhtf returns a tranformation matrix based on the modified
%Denavit-Hartenberg Parameters
%   this function returns a transformation matrix based on input of
%   modified Denavit-Hartenberg parameters as described in Dynamics of
%   Parallel robots - From Rigid Bodies to Flexible Elements (Briot and
%   Khalil) pg 40. This takes an input of a gamma, b, alpha, d, theta, and
%   r and returns a 4x4 transformation matrix.

    T = [cos(gamma)*cos(theta)-sin(gamma)*cos(alpha)*sin(theta), -cos(gamma)*sin(theta)-sin(gamma)*cos(alpha)*cos(theta),  sin(gamma)*sin(alpha), d*cos(gamma)+r*sin(gamma)*sin(alpha);
         sin(gamma)*cos(theta)+cos(gamma)*cos(alpha)*sin(theta), -sin(gamma)*sin(theta)+cos(gamma)*cos(alpha)*cos(theta), -cos(gamma)*sin(alpha), d*sin(gamma)-r*cos(gamma)*sin(alpha);
         sin(alpha)*sin(theta),                                   sin(alpha)*cos(theta),                                   cos(alpha),            r*cos(alpha)+b;
         0,                                                       0,                                                       0,                     1];
    
    for i=1:4
        for j = 1:4
            [C,c] = coeffs(T(i,j));
            for k=1:length(c)
               if C(k) < 1e-5
                   C(k) = 0;
               end
            end
            if ~isempty(c)
                T(i,j) = dot(C,c);
            end
        end
    end
end

