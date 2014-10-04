function [ x,y,w ] = maxormin( a,b,c )
%maxormin(a,b,c) finds the maximum or minimum value of the polynomial ax^2
%+ bx + c. The outputs are as follows: x is the x value of the maximum or
%minimum; y is the value of the function's maximum or minimum; w is 1 if
%there is a maximum and w is 2 if there is a minimum. 

polynom = [a b c];

% Let's find the maximum/minimum by taking the derivative and taking the
% roots

% We note that if a = 0, there is no maximum or minimum
if a == 0
    fprintf('This equation is not that of a parabola. There is no maximum or minimum.\n');
    return;
else
    x = roots(polyder(polynom));
    y = a*x^2 + b*x + c;
end

% Let's find if this is a maximum or a minimum by taking the second
% derivative and determining if it's positive or negative

if polyder(polyder(polynom)) > 0 %This corresponds to a minimum
    w = 2;
else                             %This corresonds to a maximum
    w = 1;
end

end

