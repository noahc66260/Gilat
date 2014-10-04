function y = cubic(P)
%cubic.m calculates the cube root of P by using Newton's methods.

y = P;  % first guess
n = 1;
error = 1;

while error > 0.00001
    y(n+1) = y(n) - (y(n)^3 - P)/(3*y(n)^2);
    error = abs((y(n+1) - y(n))/y(n));
    n = n + 1;
    
    if n > 3000
        fprintf('More than 3000 iterations needed.\n')
        break
    end
end

y = y(end);


end

