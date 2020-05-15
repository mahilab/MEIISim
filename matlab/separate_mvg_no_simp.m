% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/11/2016

function [M,V,G] = separate_mvg_no_simp(MVG,Qdd,g)
% =========================================================================
% Extracts the mass matrix M, the vector of centrifugal and Coriolis terms
% V, and the vector of gravity terms G from the clumped symbolic expression
% MVG. Qdd is the [n x 1] vector of joint accelerations and g is gravity.
% =========================================================================
n = length(MVG);
M = sym(zeros(n,n));
G = sym(zeros(n,1));
for i = 1:n
    for j = 1:n
        M(i,j) = diff(MVG(i),Qdd(j));
    end
    G(i,1) = diff(MVG(i),g) * g;
end
V = subs(MVG,[Qdd.',g],[zeros(size(Qdd)).',0]);
end
