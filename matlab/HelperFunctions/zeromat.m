function [ zeroed_mat ] = zeromat( mat, threshold )
%zeromat removes all coefficients that are under a specific threshold on
%   Detailed explanation goes here
    for i=1:size(mat,1)
        for j = 1:size(mat,2)
            if(~isnumeric(mat(i,j)))
                [C,c] = coeffs(expand(mat(i,j)));
                for k=1:length(C)
                   if abs(C(k)) < threshold
                       C(k) = 0;
                   end
                end
                if ~isempty(C)
                    mat(i,j) = dot(C,c);
                end
            end
        end
    end
    zeroed_mat = mat;
end