function area = triangle(a,b,c)
%triangle.m determines the area of a triangle when the lengthos of the
%sides are given
%
% Inputs:
%   a,b,c   side lengths of the triangle
%
% Outputs:
%   area    area of the triangle

s = (a + b + c)/2;

area = sqrt(s.*(s-a).*(s-b).*(s-c));

end

