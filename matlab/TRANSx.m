function [T] = TRANSx(d)

% returns Homogeneous Transform about x given d

T = [ 1 0 0 d;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];