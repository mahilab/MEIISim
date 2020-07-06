function [ T ] = mdhtf2( gamma, b, alpha, d, theta, r)
%mdhtf returns a tranformation matrix based on the modified
%Denavit-Hartenberg Parameters
%   this function returns a transformation matrix based on input of
%   modified Denavit-Hartenberg parameters as described in Dynamics of
%   Parallel robots - From Rigid Bodies to Flexible Elements (Briot and
%   Khalil) pg 40. This takes an input of a gamma, b, alpha, d, theta, and
%   r and returns a 4x4 transformation matrix.
    cosGamma = cos(gamma);
    cosTheta = cos(theta);
    cosAlpha = cos(alpha);
    sinGamma = sin(gamma);
    sinTheta = sin(theta);
    sinAlpha = sin(alpha);

    T = [cosGamma*cosTheta-sinGamma*cosAlpha*sinTheta, -cosGamma*sinTheta-sinGamma*cosAlpha*cosTheta,  sinGamma*sinAlpha, d*cosGamma+r*sinGamma*sinAlpha;
         sinGamma*cosTheta+cosGamma*cosAlpha*sinTheta, -sinGamma*sinTheta+cosGamma*cosAlpha*cosTheta, -cosGamma*sinAlpha, d*sinGamma-r*cosGamma*sinAlpha;
         sinAlpha*sinTheta,                                   sinAlpha*cosTheta,                                   cosAlpha,            r*cosAlpha+b;
         0,                                                       0,                                                       0,                     1];
     
end

