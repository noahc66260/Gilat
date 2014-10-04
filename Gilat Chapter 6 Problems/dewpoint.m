function Td = dewpoint(T,RH)
% dewpoint.m calculates the approximate dew point temperature given the
% actual temperature and relative humidity.
%
% Inputs:
%   T   actual temperature in degrees Celsius
%   RH  relative humidity in percent
%
% Outputs:
%   Td  dew point temperature in degrees Celsius

a = 17.27;
b = 237.7;  % Degrees C

f = @(T,RH) a*T./(b + T) + log(RH/100);

Td = b*f(T,RH)./(a - f(T,RH));

end

