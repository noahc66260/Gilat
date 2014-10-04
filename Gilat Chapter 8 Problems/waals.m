function V = waals(p,T,n,a,b)
%waals takes the pressure, temperature, number of moles, values of a and b
%and determines the volume of gas using Van der Waal's equation by plugging
%in these parameters and solving for V. Only appropriate values of V that
%solve the polynomial are used. 

R = 0.08206; % (L*atm)/(mole*K)
polyV = [1 -(n*b + n*R*T/p) (n^2*a/p) -(n^3*a*b/p)];
root = roots(polyV);
for i = 1:length(root)
    if imag(root(i)) == 0 && root(i) >= 0
        V = root(i);
    end
end

end

