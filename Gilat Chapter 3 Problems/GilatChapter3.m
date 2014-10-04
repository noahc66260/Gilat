% GilatChapter3.m
% Noah Ruderman
% This script file includes all the problems from Chapter 3 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 3: Mathematical Operations with Arrays
% 5.22.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
y = @(x) (2*x.^2 - 5*x + 4).^3./x.^2;
x = -2:5;

q = y(x);

for i = 1:length(x)
    fprintf('\tFor x = %+2d, y = %9.4f\n', x(i), q(i))
end
fprintf('\n\n')





%% Problem 2
disp('Problem 2')
y = @(t) 5*sqrt(t) - (t + 2).^2./(0.5*(t + 1)) + 8;
t = 0:8;

q = y(t);

for i = 1:length(t)
    fprintf('\tFor t = %d, y = %+7.4f\n', t(i), q(i))
end
fprintf('\n\n')





%% Problem 3
disp('Problem 3')
g = 9.81;   % m/s^2
h_0 = 2;    % m  initial height the ball is dropped from

h_max = @(n) h_0*0.85.^(2*n);
n = 1:8;
t = h_max(n);

for i = 1:length(t)
    fprintf('\tAfter %d bounces, the maximum height the ball reaches is %.4f m\n', ...
        n(i), t(i))
end
fprintf('\n\n')





%% Problem 4
disp('Problem 4')
g = 9.81;   % m/s^2
C_d = 0.5;
rho = 1.2;  % kg/m^3
m = 0.624;  % kg
r = 0.117;  % m
A = pi*r^2;     % m^2

t = 0:10;

v = @(t) sqrt(2*m*g/(rho*A*C_d))*(1 - exp(-sqrt(rho*g*C_d*A/(2*m))*t));
q = v(t);

for i = 1:length(t)
    fprintf('\tFor t = %2d, v = %8.4f m/s \n', t(i), q(i))
end
fprintf('\n\n')





%% Problem 5
disp('Problem 5')


disp('a.')
u = [14 25 -10];
mag_u = sqrt(u(1)^2 + u(2)^2 + u(3)^2);
fprintf('\tThe magnitude of u is %.4f \n', mag_u)


disp('b.')
mag_u = sqrt(sum(u.^2));
fprintf('\tThe magnitude of u is %.4f \n', mag_u)
fprintf('\n\n')





%% Problem 6
disp('Problem 6')
g = 9.81;   % m/s^2
v_0 = 100;  % m/s
theta = 79; % degrees

t = 0:2:20; % s

x = @(theta,t) v_0*cosd(theta).*t;
y = @(theta,t) v_0*sind(theta).*t - 1/2*g*t.^2;
r = @(theta,t) sqrt(x(theta,t).^2 + y(theta,t).^2);

q = r(theta,t);

for i = 1:length(q)
    fprintf('\tFor t = %2d s, r = %9.4f m \n', t(i), q(i))
end
fprintf('\n\n')





%% Problem 7
disp('Problem 7')


disp('a.')
u = [4 9 -5];
v = [-3 6 -7]';

u*v


disp('b.')
u = [4 9 -5];
v = [-3 6 -7];

dot(u,v)





%% Problem 8
disp('Problem 8')
x = 2:2:10;
y = 3:3:15;

z = (y./x).^2 + (x + y).^((y - x)./x)





%% Problem 9
disp('Problem 9')
h = 0.7;
k = 8.85;
x = 1:5;
y = 2.1:-.1:1.7;
z = 2:.5:4;

G = (h*x + k*y)./(x + y).^h + exp(h*y./z)./z.^(y./x)





%% Problem 10
disp('Problem 10')
format long
x = [1 0.5 0.1 0.01 0.001 0.00001 0.0000001];
q = (exp(x) - 1)./x;

for i = 1:length(x)
    fprintf('\tFor x = %.7f, (exp(x) - 1)/x = %.9f\n', x(i), q(i))
end
fprintf('\n\n')
format short





%% Problem 11
disp('Problem 11')
fprintf('\tFor reference, the value of pi is %.4f\n', pi)

disp('a.')
n = 0:100;
q = (-1).^n ./(2*n + 1);
t = 4*sum(q);
fprintf('\tFor n = %7d, the partial sum is %.4f\n', n(end), t)


disp('b.')
n = 0:10000;
q = (-1).^n ./(2*n + 1);
t = 4*sum(q);
fprintf('\tFor n = %7d, the partial sum is %.4f\n', n(end), t)

disp('c.')
n = 0:1000000;
q = (-1).^n ./(2*n + 1);
t = 4*sum(q);
fprintf('\tFor n = %7d, the partial sum is %.4f\n', n(end), t)
fprintf('\n\n')





%% Problem 12
disp('Problem 12')
fprintf('\tFor reference, the value of ln2 is %.4f\n', log(2))


disp('a.')
n = 0:50;
q = 1./((2*n+1).*(2*n+2));
t = sum(q);
fprintf('\tFor n = %4d, the partial sum is %.4f\n', n(end), t)


disp('b.')
n = 0:500;
q = 1./((2*n+1).*(2*n+2));
t = sum(q);
fprintf('\tFor n = %4d, the partial sum is %.4f\n', n(end), t)


disp('c.')
n = 0:5000;
q = 1./((2*n+1).*(2*n+2));
t = sum(q);
fprintf('\tFor n = %4d, the partial sum is %.4f\n', n(end), t)
fprintf('\n\n')





%% Problem 13
disp('Problem 13')
L_max = 50;     % cm
tau = 0.5;  % years
K = [0.25 .5 .75];  % yrs^(-1)
t = 2;

L = L_max*(1 - exp(-K*(t + tau)));

for i = 1:length(K)
    fprintf('\tAt 2 years of age, for K = %.2f, the length of a fish is %.4f cm\n', ...
        K(i), L(i))
end
fprintf('\n\n')





%% Problem 14
disp('Problem 14')
A = [ 5 2 4; 1 7 -3; 6 -10 0];
B = [ 11 5 -3; 0 -12 4; 2 6 1];
C = [7 17 1; 10 3 -2; 8 -5 9];

disp('a.')
A + B
B + A
fprintf('\tWe see that addition of matrices is commutative\n')


disp('b.')
A + (B + C)
(A + B) + C
fprintf('\tWe see that addition of matrices is associative\n')


disp('c.')
5*(A + C)
5*A + 5*C
fprintf('\tWe see that multiplication of matrices by scalars is distributive\n')


disp('d.')
A*(B + C) 
A*B + A*C
fprintf('\tWe see that matrix multiplication is distributive\n')
fprintf('\n\n')





%% Problem 15
disp('Problem 15')
A = [ 5 2 4; 1 7 -3; 6 -10 0];
B = [ 11 5 -3; 0 -12 4; 2 6 1];
C = [7 17 1; 10 3 -2; 8 -5 9];


disp('a.')
if (all(all( A*B == B*A )))
    fprintf('\tA*B == B*A\n')
else
    fprintf('\tA*B ~= B*A\n')
end

    
disp('b.')
if (all(all( A*(B*C) == (A*B)*C )))
    fprintf('\tA*(B*C) == (A*B)*C \n')
else
    fprintf('\tA*(B*C) ~= (A*B)*C \n')
end

disp('c.')
if (all(all( (A*B)' == B'*A' )))
    fprintf('\t(A*B)'' == B''*A'' \n')
else
    fprintf('\t(A*B)'' ~= B''*A'' \n')
end


disp('d.')
if (all(all( (A + B)' == A' + B' )))
    fprintf('\t(A + B)'' == A'' + B'' \n')
else
    fprintf('\t(A + B)'' ~= A'' + B'' \n')
end
fprintf('\n\n')





%% Problem 16
disp('Problem 16')
v1 = 680;   % m/s
theta1 = 65;    % degrees
v2 = 780;   % m/s
theta2 = 42;    % degrees
g = 9.81;   % m/s^2

v_x1 = v1*cosd(theta1);
v_y1 = v1*sind(theta1);
v_x2 = v2*cosd(theta2);
v_y2 = v2*sind(theta2);

if v_y1 < v_y2
    fprintf('\tProjectile A will hit the ground first\n')
elseif v_y1 > v_y2
    fprintf('\tProjectile B will hit the ground first\n')
else
    fprintf('\tBoth projectiles will hit the ground at the same time\n')
end

t_f = 2*v_y2/g;
t = linspace(0, t_f, 11);

x_distance = @(t) t*(v_x1 - v_x2);
y_distance = @(t) t*(v_y1 - v_y2);
r_distance = @(t) sqrt(x_distance(t).^2 + y_distance(t).^2);

q = r_distance(t);

for i = 1:length(t)
    fprintf('\tFor t = %8.4f s, the distance between A and B is %10.4f m\n', ...
        t(i), q(i))
end
fprintf('\n\n')





%% Problem 17
disp('Problem 17')


% disp('a.')
u = 0:0.05:1;


% disp('b.')
k = 0.25;
p = k*u.*(1 - u)./(k + u);


% disp('c.')
max_p05 = max(p);


disp('d.')
u = 0:0.01:1;
k = 0.25;
p = k*u.*(1 - u)./(k + u);
max_p01 = max(p);

error = abs((max_p01 - max_p05)/max_p05)*100;
fprintf('\tThe relative error is %.4f \n\n\n', error)





%% Problem 18
disp('Problem 18')
A = [1.5, -2, 1, 3, 0.5;
    3, 1, -1, 4, -3;
    2, 6, -3, -1, 3;
    5, 2, 4, -2, 6;
    -3, 3, 2, 5, 4];
b = [7.5, 16, 78, 71, 54]';

x = A\b;

fprintf('\tx = %+6.4f\n', x(1))
fprintf('\ty = %+6.4f\n', x(2))
fprintf('\tz = %+6.4f\n', x(3))
fprintf('\tu = %+6.4f\n', x(4))
fprintf('\tw = %+6.4f\n', x(5))
fprintf('\n\n')





%% Problem 19
disp('Problem 19')
V1 = 38; V2 = 20; V3 = 24;    % V
R1 = 15; R2 = 18; R3 = 10; R4 = 9; R5 = 5; R6 = 14; R7 = 8; R8 = 13;    % ohms

A = [-(R1 + R2)         (R2)                  0              0              0       ;
       (R2)      -(R2 + R3 + R4 + R7)        (R3)           (R4)           (R7)     ;
        0               (R3)           -(R3 + R5 + R6)      (R5)           (R6)     ;
        0               (R4)                 (R5)        -(R4 + R5)         0       ;
        0               (R7)                 (R6)            0       -(R6 + R7 + R8)];

b = [-V1, -V2, 0, V3, V1]';

x = A\b;

I = zeros(1,8);
I(1) = x(1);
I(2) = x(1) - x(2);
I(3) = x(2) - x(3);
I(4) = x(2) - x(4);
I(5) = x(4) - x(3);
I(6) = x(3) - x(5);
I(7) = x(5) - x(2);
I(8) = x(5);

I = abs(I);

for i = 1:length(I)
    fprintf('\tThe current through resistor %d is %6.4f amps\n', i, I(i))
end
fprintf('\n\n')




