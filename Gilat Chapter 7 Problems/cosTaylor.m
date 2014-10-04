function y = cosTaylor(x)
%cosTaylor.m calculates the cosine of an angle in degrees by summing terms
%from the Taylor series until the error is less than 0.000001.

x = x*pi/180;
newsum = 0;
oldsum = 1;
error = 1;
n = 1;
while error > 0.000001
    newsum = oldsum + (-1)^(n)/fact(2*n)*x^(2*n);
    error = abs((newsum - oldsum)/oldsum);
    oldsum = newsum;
    n = n + 1;
end

y = newsum;



function y = fact(x)
%fact.m calculates the factorial of the scalar x. 

if x < 0 || int32(x) ~= x
    fprintf('Error, input must be a nonpositive integer.\n')
    return
    y = -1;
else
    if x ~= 0
        y = x*fact(x-1);
    else
        y = 1;
    end
end


end


end

