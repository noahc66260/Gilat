function dxyzdt = dladderdt(t, q)
%dladderdt.m will be used to solve the ODE for calculating the position of
%the ladder in question 3 of Chapter 9 in Gilat's book. 

x = q(1);
y = q(2); 
z = q(3);


dxdt = (0.6)*x./sqrt(x.^2 + y.^2 + z.^2) ...
    - (pi/36)*z.*x./sqrt(x.^2 + y.^2) ...
    - (2*pi/45)*y;

dydt = (0.6)*y./sqrt(x.^2 + y.^2 + z.^2) ...
    - (pi/36)*z.*y./sqrt(x.^2 + y.^2) ...
    + (2*pi/45)*x;

dzdt = (0.6)*z./sqrt(x.^2 + y.^2 + z.^2) ...
    + (pi/36)*sqrt(x.^2 + y.^2);

dxyzdt(1,1) = dxdt;
dxyzdt(2,1) = dydt;
dxyzdt(3,1) = dzdt;


end

