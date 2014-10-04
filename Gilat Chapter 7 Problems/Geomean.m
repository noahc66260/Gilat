function GM = Geomean( x )
%Geomean.m calculates the geometric mean of the elements in a vector.

GM = 1;
n = length(x);
for i = 1:n
    GM = GM * x(i);
end
GM = GM^(1/n);


end

