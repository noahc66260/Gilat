function RV = bandpass( R,C,L,w )
%bandpass.m calculates the ratio of voltages (V_0 / V_i) in a bandpass
%filter.
%   Input:
%       R = the size of the reisistor in ohms
%       C = the size of the capacitor in farads
%       L = the inductance of the coil in Henrys
%       w= the frequency of the input signal in radians/second
%   Output:
%       RV = the ratio of the voltages

RV = w.*R.*C ./ sqrt((1-w.^2.*L.*C).^2 + (w.*R.*C).^2);


end

