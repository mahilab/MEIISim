function [T] = TRANSy(d)

% returns Homogeneous Transform about y given d

T = [ 1 0 0 0;
      0 1 0 d;
      0 0 1 0;
      0 0 0 1];