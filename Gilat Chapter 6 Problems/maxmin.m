function [ x,y ] = maxmin( a,b,c )
%maxormin.m calculates that maximum or minimum of a quadratic equation of
%the form y = ax^2 + bx + c.
%   a = coefficient of the x^2 term in the quadratic equation
%   b = coefficient of the x^1 term in the quadratic equation
%   c = coefficient of the x^0 term in the quadratic equation
%   x = x coordinate of the max or min
%   y = function value at the max or min
% We determine the max or min by where the derivative of the function is 0

x = -b/(2*a);
quad_eq = [a b c];
y = polyval(quad_eq,x);


end

