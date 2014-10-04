function [Smax, Smin] = princstress(Sxx, Syy, Sxy)
% princstress.m finds the principal stresses for a two-dimensional state of
% stress defined by problem 17 of Gilat's MATLAB: An Introduction with
% Applications in Chapter 6.
% 
% Inputs:
%   Sxx     
%   Syy     
%   Sxy     
%
% Outputs:
%   Smax    maximum stress
%   Smin    minimum stress

Smax = (Sxx + Syy)/2 + sqrt(((Sxx - Syy)/2).^2 + Sxy.^2);
Smin = (Sxx + Syy)/2 - sqrt(((Sxx - Syy)/2).^2 + Sxy.^2);

end

