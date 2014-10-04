function M = amort(P,r,N)
%amort.m calculates the monthly payment of a loan
% Inputs:
%   P   loan amount
%   r   annual interest rate in percent
%   N   length of the loan in years
%
% Outputs:
%   M   monthly payment

M = P*(r/1200)/(1 - (1 + r/1200)^(-12*N));

end

