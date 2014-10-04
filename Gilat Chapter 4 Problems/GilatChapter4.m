% GilatChapter4.m
% Noah Ruderman
% This script file includes all the problems from Chapter 4 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 4: Using Script Files and Managing Data
% 5.15.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
% So we know that R2 = 1.25*R1, V = 250 cm^3, 
% V = 1/3*pi*h*(R1^2 + R2^2 + R1*R2)
% S = pi*(R1 + R2)*sqrt((R2 - R1)^2 + h^2) + pi*R1^2
% We can write R1 as a function of V,h

V = 250;     % cm^3
h = 5:10;   % cm
R1 = sqrt(3*V./(pi*h*(2.25 + 1.25^2)));
R2 = 1.25*R1;
S = pi*(R1 + R2).*sqrt((R2 - R1).^2 + h.^2) + pi*R1.^2;

for i = 1:length(h)
    fprintf('\tFor a height of %2d cm, R1 = %.4f cm, R2 = %.4f cm, ',...
        h(i), R1(i), R2(i));
    fprintf('and S = %.4f cm^2.\n', S(i));
end
fprintf('\n\n')





%% Problem 2
disp('Problem 2')
% I will calculate the angle by subtracting the angle of the viewer to the
% lower end of the screen from the horizontal from the angle of the viewer
% to the upper end of the screen from the horizontal

x = [30 45 60 75 90];   % ft
viewer_height = x*tand(8);  % ft
screen_height_lower_end = 6;    % ft
screen_height_upper_end = 6 + 24;   % ft
theta_lower = atand((screen_height_lower_end - viewer_height)./x);  % degrees
theta_upper = atand((screen_height_upper_end - viewer_height)./x);  % degrees
theta = theta_upper - theta_lower; % degrees

for i = 1: length(x)
    fprintf('\tFor a viewer sitting %2d feet from the screen, ', x(i))
    fprintf('theta = %6.4f degrees\n', theta(i))
end
fprintf('\n\n')





%% Problem 3
disp('Problem 3')
v_run = 3; % m/s
v_swim = 1; % m/s
L = 48; % m
d_s = 30; % m
d_w = 42; % m
y = 1:48; % m

% we know that time = distance/velocity

d1 = sqrt(y.^2 + d_s^2);
d2 = sqrt((L-y).^2 + d_w^2);
t = d1./v_run + d2./v_swim;
[t_min index] = min(t);
y_min = y(index);

fprintf('\tThe minimum time is %.4f seconds and the entry point y is ', t_min)
fprintf('%.0f meters.\n', y_min)

fprintf('\tsin(phi)/sin(alpha) = %.4f\n', sin(y_min/d1(index))/sin((L-y_min)/d2(index)))
fprintf('\tv_run/v_swim        = %.4f\n', v_run/v_swim)
fprintf('\tWe see that Snell''s Law is not satisfied in this situation.\n')
fprintf('\n\n')





%% Problem 4
disp('Problem 4')
t_halflife = 6;     % hrs
k = -log(2)/t_halflife; 

t = 0:2:24;     % hrs
relative_Tc99 = exp(k*t);

for i = 1:length(t)
    fprintf('\tAfter %2d hours, the relative amount of technetium-99 left is', t(i))
    fprintf(' %.4f\n', relative_Tc99(i));
end
fprintf('\n\n')




%% Problem 5
disp('Problem 5')
n = 1:10;   % years
A = 1000;   % USD initial investment
r = 6.5;   % annual compound interest rate in percent

B = A*(1 + r/100).^n;

format bank
table(:,1) = n';
table(:,2) = B';
disp('          Year       Balance')
disp(' ')
disp(table)
fprintf('\n\n')
format short





%% Problem 6
disp('Problem 6')
t = 1:10;   % seconds
a = 1.55;   % m/s^2 acceleration

v = a*t;    % m/s
d = 1/2*a*t.^2;     % m

fprintf(' Time (s)  Distance (m)  Velocity (m/s) \n\n')
for i = 1:length(t)
    fprintf(' %4d      %8.4f      %8.4f \n', t(i), d(i), v(i))
end
fprintf('\n\n')





%% Problem 7
disp('Problem 7')
clear
a = 34172; % I think this may be a typo in the book
b = 7.9622; 

T_C = 0:2:42;
T_K = T_C + 273;

p = 10.^(b - 0.05223*a./T_K);

fprintf(' Temperature (Celsius)      Pressure (mmHg) \n')
for i = 1:length(p)
    fprintf('         %-10d   %16.4f   \n', T_C(i), p(i))
end
fprintf('\n\n')





%% Problem 8
disp('Problem 8')
T = 200:20:400;
C_p_SO2 = polyval([8.606*10^(-9), -3.105*10^(-5), 3.904*10^(-2), 38.91], T);
C_p_SO3 = polyval([32.40*10^(-9), -8.540*10^(-5), 9.188*10^(-2), 48.50], T);
C_p_O2  = polyval([1.311*10^(-9), -0.6076*10^(-5), 1.158*10^(-2), 29.10], T);
C_p_N2  = polyval([-2.871*10^(-9), -0.5723*10^(-5), 0.2199*10^(-2), 29.00], T);

fprintf(' Temperature (Celsius)           Heat Capacity (J/(g mol C)) \n')
fprintf('                              SO_2     SO_3      O_2      N_2     \n')
for i = 1:length(T)
    fprintf('   %10d                %6.4f   %6.4f   %6.4f   %6.4f \n', ...
        T(i), C_p_SO2(i), C_p_SO3(i), C_p_O2(i), C_p_N2(i))
end
fprintf('\n\n')





%% Problem 9
disp('Problem 9')
% We know that Ax = b where A is a matrix of the coefficients of the linear
% system of equations and b is a vector with the solutions. We find x by
% finding A^(-1)*b

T = [25 150 300];
C_p_SO2 = polyval([8.606*10^(-9), -3.105*10^(-5), 3.904*10^(-2), 38.91], T);
C_p_SO3 = polyval([32.40*10^(-9), -8.540*10^(-5), 9.188*10^(-2), 48.50], T);
C_p_O2  = polyval([1.311*10^(-9), -0.6076*10^(-5), 1.158*10^(-2), 29.10], T);
C_p_N2  = polyval([-2.871*10^(-9), -0.5723*10^(-5), 0.2199*10^(-2), 29.00], T);
b = [39.82 44.72 49.10 1]';

A = zeros(4,4);
A(:,1) = [C_p_SO2 1]';
A(:,2) = [C_p_SO3 1]';
A(:,3) = [C_p_O2 1]';
A(:,4) = [C_p_N2 1]';

x = A\b;

fprintf('\tThe fraction of SO_3 is %.4f \n', x(1))
fprintf('\tThe fraction of SO_2 is %.4f \n', x(2))
fprintf('\tThe fraction of O_2 is %.4f \n', x(3))
fprintf('\tThe fraction of N_2 is %.4f \n', x(4))
fprintf('\n\n')





%% Problem 10
disp('Problem 10')

v_s = input('Enter the source voltage: ');
R = input('Enter the values of the resistors in vector form (ex. [12 34 93]): \n');
i = v_s./R;
P = v_s*i;

fprintf('  Resistance (Ohms)     Current     Power \n')
for j = 1:length(R)
    fprintf('  %-17d     %-7.4f     %-.4f   \n', R(j), i(j), P(j));
end

i_s = v_s*sum(1./R);
P_total = sum(P);

fprintf('The source current is %.4f amps.\n', i_s)
fprintf('The total power is %.4f watts.\n', P_total)
fprintf('\n\n')





%% Problem 11
disp('Problem 11')

A = [cosd(45)       1           0           0           0           0           0           0           0;
     -cosd(45)      0           0           1       cosd(48.81)     0           0           0           0;
     -sind(45)      0           -1          0       -sind(48.81)    0           0           0           0;
     0              0           0           -1          0           0           0         cosd(48.81)   0;
     0              0           0           0           0           0           -1       -sind(48.81)   0;
     0              0           0           0       -cosd(48.81)    -1          0           0           1;
     0              0           0           0        sind(48.81)    0           1           0           0;
     0              0           0           0           0           0           0         sind(48.81)   0;
     0              0           0           0           0           0           0        -cosd(48.81)  -1;];
 
b = [0; 0; 1000; 0; 500; 0; 4000; -1107.14; 0];

x = A\b;

fprintf(' Force   Magnitude (N) \n')
for i = 1:length(x)
    fprintf('   %-2d      %-+9.4f  \n', i, x(i))
end
fprintf('\n\n')





%% Problem 12
disp('Problem 12')
x1 = -3;
y1 = 6.8;

x2 = -1.5;
y2 = 15.2;

x3 = 0.5;
y3 = 14.5;

x4 = 2;
y4 = -21.2;

x5 = 5;
y5 = 10;

A = [x1^4     x1^3      x1^2      x1        1;
     x2^4     x2^3      x2^2      x2        1;
     x3^4     x3^3      x3^2      x3        1;
     x4^4     x4^3      x4^2      x4        1;
     x5^4     x5^3      x5^2      x5        1;];

b = [y1; y2; y3; y4; y5];

x = A\b;

fprintf('\ta = %.4f \n', x(1))
fprintf('\tb = %.4f \n', x(2))
fprintf('\tc = %.4f \n', x(3))
fprintf('\td = %.4f \n', x(4))
fprintf('\te = %.4f \n', x(5))
fprintf('\n\n')





%% Problem 13
disp('Problem 13')
A = [1  5   2   0;
     2  3   1   1;
     0  3   3   2;
     1  4   2   1];

b = [18; 15; 0; 12];

x = A\b;
x = int32(x);

fprintf('\tAn Eagle is worth %d points.\n', x(1))
fprintf('\tA Birdie is worth %d points.\n', x(2))
fprintf('\tA Bogey subtracts %d points.\n', -x(3))
fprintf('\tA double bogey subtracts %d points.\n', -x(4))
fprintf('\n\n')





%% Problem 14
disp('Problem 14')
A = [1  0   0   -1   0   0   0;
     1  0   0   0   -1   0   0;
     0  1   0   0    0  -1   0;
     0  3   0   0   -4  -1  -1;
     0  0   1   0    0   0  -2;
     0  -1   1  -2   2   0   0;
     1  0   0   0    0   0   0];
 
b = [0; 0; 0; 0; 0; 0; 3];

x = A\b;
x = int32(x);

fprintf('\ta = %d\n', x(1))
fprintf('\tb = %d\n', x(2))
fprintf('\tc = %d\n', x(3))
fprintf('\td = %d\n', x(4))
fprintf('\te = %d\n', x(5))
fprintf('\tf = %d\n', x(6))
fprintf('\tg = %d\n', x(7))
fprintf('\n\n')




