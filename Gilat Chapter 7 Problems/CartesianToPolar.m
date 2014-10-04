function [theta radius] = CartesianToPolar(x,y)
%CartesianToPolar.m converts cartesian to polar coordinates where the angle
%theta is in degrees.

radius = sqrt(x.^2 + y.^2);
theta = atand(y/x);


end

