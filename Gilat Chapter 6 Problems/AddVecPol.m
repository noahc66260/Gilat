function [r th] = AddVecPol(r1, th1, r2, th2)
%AddVecPol.m adds two vectors whose components are given in polar
%coordinates.
%
% Inputs:
%   r1  length of first vector
%   th1     angle of first vector in degrees
%   r2  length of second vector
%   th2     angle of second vector in degrees 
%
% Outputs:
%   r   length of new vector
%   th  angle of new vector in degrees

x1 = r1.*cosd(th1);
y1 = r1.*sind(th1);
x2 = r2.*cosd(th2);
y2 = r2.*sind(th2);

x = x1 + x2; 
y = y1 + y2;

r = sqrt(x.^2 + y.^2);
th = atand(y./x);

end

