function [ cpm ] = crossprodmat( vector )
%crossprodmat takes in a vector, and returns the cross-product matrix
%associated with the vector
%   input is a 3-vector [rx, ry, rz], and the output is a 3x3 matrix form 
%   as follows:
%             [0  -rz  ry;
%              rz  0  -rx;
%             -ry  rx  0 ]

rx = vector(1);
ry = vector(2);
rz = vector(3);

cpm = [0  -rz  ry;
       rz  0  -rx;
      -ry  rx  0 ];

end

