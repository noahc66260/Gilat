function pr = prime(n)
%prime.m finds all prime numbers between 1 and n and creates a vector with
%these values. An error message is given if the input is negative or not an
%integer.

if n < 0 || int32(n) ~= n
    fprintf('The input argument must be a positive integer.\n\n')
    return
end

pr = [];
for i = 2:n
    isprime = 1;
    for j = 2:floor(i/2)
        if rem(i,j) == 0
            isprime = 0;
        end
    end
    
    if isprime
        pr = [pr i];
    end
end


end

