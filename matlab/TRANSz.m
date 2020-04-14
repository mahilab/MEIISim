function [T] = TRANSz(d)

% returns Homogeneous Transform about z given d

T = [ 1 0 0 0;
      0 1 0 0;
      0 0 1 d;
      0 0 0 1];