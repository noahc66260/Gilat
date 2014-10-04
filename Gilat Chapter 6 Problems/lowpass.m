function RV = lowpass(R,C,w)
% lowpass.m calculates the ratio of the magnitude of the voltages
%
% Inputs:
%   R   the size of the resistor (ohms)
%   C   the size of the capacitor in Farads
%   w   the frequency of the input signal in rad/s
%
% Outputs:
%   RV  the ratio of the magnitude of the voltages

RV = 1./sqrt(1 + (w.*R.*C).^2);

end

