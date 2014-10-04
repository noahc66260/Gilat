function fx = functionFx( x )
%functionFx stores the values of the image of a function when a subset of
%its domain is specified as the input argument.

fx = zeros(size(x));
for i = 1:length(x)
    if x(i) <= -1
        fx(i) = 15;
    elseif x(i) <= 1
        fx(i) = -5*x(i) + 10;
    elseif x(i) <= 3
        fx(i) = -10*x(i)^2 + 35*x(i) - 20; 
    elseif x(i) <= 4
        fx(i) = -5*x(i) + 10;
    else
        fx(i) = -10;
    end
end


end

