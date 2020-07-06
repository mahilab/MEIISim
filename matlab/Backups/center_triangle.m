function [ new_point ] = center_triangle( point1, point2, point3)
%center_trangle Summary of this function goes here
%   Detailed explanation goes here
%     syms t1 t2
%     mid1 = (point1 + point2)/2;
%     mid2 = (point2 + point3)/2;
%     line1 = mid1 + (point3 - mid1)*t1
%     line2 = mid2 + (point1 - mid2)*t2
%     [T1, ~] = solve(line1 == line2, t1, t2);
%     new_point = simplify(expand(subs(line1, T1)));
    line1 = mid1 + (point3 - mid1)*t1
    
end
    