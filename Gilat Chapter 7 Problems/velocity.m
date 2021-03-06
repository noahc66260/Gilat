function v = velocity(t)
%velocity.m calculates the velocity of a particle that obeys the piecewise
%position time dependent functions outlined in problem 23, chapter 7 in
%MATLAB: An Introduction with Applications by Amos Gilat

x1 = [     0      0  0.5    0      0];     % for 0 <= t <= 10 s
x2 = [     0   0.05   -1   15    -50];  % for 10 <= t <= 20 s
x3 = [0.0025  -0.15    0  135  -1650]; % for 20 <= t <= 30 s

v1 = polyder(x1);
v2 = polyder(x2);
v3 = polyder(x3);

if t >= 0 & t <= 10
    v = polyval(v1,t);
elseif t <= 20
    v = polyval(v2,t);
elseif t <= 30
    v = polyval(v3,t);
end

end

