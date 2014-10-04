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

