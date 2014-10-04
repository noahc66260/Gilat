% GilatChapter5.m
% Noah Ruderman
% This script file includes all the problems from Chapter 5 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 5: Two-Dimensional Plots
% 5.16.11

close all
clear
clc
format loose

%% Problem 1
% disp('Problem 1')
x1 = -4:.1:4;
x2 = -8:.1:8;
fx = [0.01 0 -0.45 0.5 -2];

figure()
subplot(1,2,1)
plot(x1,polyval(fx,x1))
xlabel('x')
ylabel('f(x)')

subplot(1,2,2)
plot(x2,polyval(fx,x2))
xlabel('x')
ylabel('f(x)')





%% Problem 2
% disp('Problem 2')
x = -10:.1:10;
fx = 5./(1+exp(5.5 - 1.5*x)) - x.^2/20;

figure()
plot(x,fx)
title('Question 2 Graph')





%% Problem 3
% disp('Problem 3')

figure()
fplot(@(x) 40/(1+(x-4).^2) + 5*sin(20*x/pi), [-10 10])
title('Question 3 Graph')





%% Problem 4
% disp('Problem 4')
x1 = -4:.05:1.7;
x2 = 2.3:.05:8;

fx = @(x) (x.^2 - 4*x - 5)./(x - 2);

y1 = fx(x1);
y2 = fx(x2);

figure()
plot(x1,y1,x2,y2)
title('Question 4 Graph')





%% Problem 5
% disp('Problem 5')
x1 = -10:.005:-2.01;
x2 = -1.99:.005:4.99;
x3 = 5.01:.005:10;

fx = @(x) (4*x - 30)./(x.^2 - 3*x - 10);

y1 = fx(x1);
y2 = fx(x2); 
y3 = fx(x3);

figure()
plot(x1,y1,x2,y2,x3,y3)
title('Question 5 Graph')
axis([-10 10 -20 20])





%% Problem 6
% disp('Problem 6')

figure()
fplot(@(x) 3*x*cos(x)^2 - 2*x, [-2*pi 2*pi])
hold on
fplot(@(x) 3*cos(x)^2 - 3*x*2*cos(x)*sin(x) - 2, [-2*pi 2*pi], ...
    'LineStyle', '--')
title('Question 6 Graph')
xlabel('x axis')
ylabel('f(x) and f''(x)')
legend('f(x)', 'f''(x)')
hold off





%% Problem 7
% disp('Problem 7')
v_s = 12;   % volts
r_s = 2.5;  % ohms

figure()
fplot(@(R_L) v_s^2*R_L/(R_L + r_s)^2, [1 10])
title('Question 7 Graph')
xlabel('Resistance (Ohms)')
ylabel('Power (Watts)')





%% Problem 8
disp('Problem 8')
t = 0:.1:4;
v_a = zeros(length(t),2);
v_a(:,1) = (14 - 6*t');
v_b = zeros(length(t),2);
v_b(:,1) = (-25*cosd(30) + 14*cosd(30)*t');
v_b(:,2) = (-25*sind(30) + 14*sind(30)*t');

distance = sqrt(dot(v_a-v_b, v_a-v_b,2)');

figure()
plot(t+7,distance)
title('Problem 8 Graph')
xlabel('Time (AM)')
ylabel('Distance (miles)')

fprintf('\tFrom the graph, people from both ships should be able to see ')
fprintf('each other between \n8:45 and 9:15 AM.\n')
fprintf('\n\n')





%% Problem 9
% disp('Problem 9')
figure()
fplot(@(x) 693.8 - 68.8*cosh(x/99.7),[-300 300])
title('Problem 9 Graph')





%% Problem 10
% disp('Problem 10')
t = 0:.1:30; % seconds
x = polyval([-0.28 6.5 61],t);
y = polyval([0.18 -8.5 65],t);

figure()
subplot(2,2,1)
plot(x,y)
title('Problem 10 Part A Graph')
xlabel('x axis')
ylabel('y axis')

xy_vec = [x' y'];
r_magnitude = sqrt(dot(xy_vec, xy_vec, 2))';

subplot(2,2,2)
plot(t,r_magnitude)
xlabel('Time (s)')
ylabel('Position Vector Magnitude (m)')
title('Part B Graph')

theta = atand(y./x)';

subplot(2,2,3)
plot(t,theta)
xlabel('Time (s)')
ylabel('Theta (degrees)')
title('Part C Graph')





%% Problem 11
% disp('Problem 11')
Table = struct();
Table.Sun =             struct('Temp',  5840, 'rel_L',      1, 'rel_R',    1);
Table.Spica =           struct('Temp', 22400, 'rel_L',  13400, 'rel_R',  7.8);
Table.Regulus =         struct('Temp', 13260, 'rel_L',    150, 'rel_R',  3.5);
Table.Alioth =          struct('Temp',  9400, 'rel_L',    108, 'rel_R',  3.7);
Table.Barnard_Star =    struct('Temp',  3130, 'rel_L', 0.0004, 'rel_R', 0.18);
Table.Epsilon_Indi =    struct('Temp',  4280, 'rel_L',   0.15, 'rel_R', 0.76);
Table.Beta_Crucis =     struct('Temp', 28200, 'rel_L',  34000, 'rel_R',    8);

% Our values for rel_L from the table
y1 = [Table.Sun.rel_L, Table.Spica.rel_L, Table.Regulus.rel_L, Table.Alioth.rel_L, ...
    Table.Barnard_Star.rel_L, Table.Epsilon_Indi.rel_L, Table.Beta_Crucis.rel_L];

% Our calculated values of rel_L using Temp and rel_R from the table
y2 = [  Table.Sun.rel_R^2*(Table.Sun.Temp/Table.Sun.Temp)^4,
        Table.Spica.rel_R^2*(Table.Spica.Temp/Table.Sun.Temp)^4,
        Table.Regulus.rel_R^2*(Table.Regulus.Temp/Table.Sun.Temp)^4,
        Table.Alioth.rel_R^2*(Table.Alioth.Temp/Table.Sun.Temp)^4,
        Table.Barnard_Star.rel_R^2*(Table.Barnard_Star.Temp/Table.Sun.Temp)^4
        Table.Epsilon_Indi.rel_R^2*(Table.Epsilon_Indi.Temp/Table.Sun.Temp)^4,
        Table.Beta_Crucis.rel_R^2*(Table.Beta_Crucis.Temp/Table.Sun.Temp)^4];

x = [Table.Sun.Temp, Table.Spica.Temp, Table.Regulus.Temp, Table.Alioth.Temp, ...
    Table.Barnard_Star.Temp, Table.Epsilon_Indi.Temp, Table.Beta_Crucis.Temp];

figure()

loglog(x,y1,'*')
hold on
loglog(x,y2,'o')
title('Problem 11 Graph')
xlabel('Log Temperature')
ylabel('Relative Luminosity')
set(gca,'XDir', 'reverse')
axis([3000 30000  .0001 100000])
legend('Experiment', 'Theory')





%% Problem 12
% disp('Problem 12')
x = [-0.1 0.8 0 10 -70];
v = polyder(x);
a = polyder(v);

t = 0:.01:8;

figure()
subplot(2,2,1)
plot(t,polyval(x,t))
xlabel('Time (s)')
ylabel('Distance (m)')
title('Problem 12 Graph Position vs. time')

subplot(2,2,2)
plot(t,polyval(v,t))
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs. time')

subplot(2,2,3)
plot(t,polyval(a,t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Acceleration vs. time')





%% Problem 13
% disp('Problem 13')
A_0 = (0.0064)^2*pi; % meters^2
L_0 = 25;   % mm

F = [0 13345 26689 40479 42703 43592 44482 44927 45372 46276 47908];
F = [F 49035 50265 53213 56161]; 

L = [25 25.037 25.073 25.113 25.122 25.125 25.132 25.144 25.164 25.208];
L = [L 25.409 25.646 26.084 27.398 29.150];

sigma_e = F/A_0;
epsilon_e = (L-L_0)/L_0;

sigma_t = F/A_0.*L/L_0;
epsilon_t = log(L/L_0);

figure()
plot(epsilon_e,sigma_e,epsilon_t,sigma_t)
title('Problem 13 Graph')
xlabel('Stress')
ylabel('Strain (Pa)')
legend('Engineering', 'True')




%% Problem 14
% disp('Problem 14')
Q1 = 4;     % Liters/min
Q2 = 5;     % Liters/min
PG = 2:60;  % mmHg
Av1 = Q1./sqrt(PG); % cm^2
Av2 = Q2./sqrt(PG); % cm^2

figure()
plot(PG,Av1,PG,Av2)
xlabel('PG (mmHg)')
ylabel('Area (cm^2)')
title('Problem 14 Graph')
legend('Q = 4 L/min', 'Q = 5 L/min')
axis([2 60 0 4])





%% Problem 15
% disp('Problem 15')
R = 4;  % Ohms
L = 1.3;    % Henrys
V = 12;     % Volts

t1 = 0:0.01:0.5;
t2 = 0.51:0.01:2;
t = [t1 t2];

i = [(V / R * (1 - exp(-R*t1/L)))   (exp(-R*t2/L) * V / R * (exp(R/2/L) - 1))];

figure()
plot(t, i)
xlabel('Time (s)')
ylabel('Current (A)')
title('Problem 15 Graph')





%% Problem 16
% disp('Problem 16')
G_inf = 5;  % ksi
c = 0.05; 
tau1 = 0.05;    % seconds
tau2 = 500;     % seconds
omega = 10.^(linspace(-4,3,100));   % seconds^(-1)

G_prime = G_inf*(1 + c/2*log((1+(omega*tau2).^2)./(1+(omega*tau1).^2)));
G_doubleprime = c*G_inf*(atan(omega*tau2)-atan(omega*tau1));

figure()
subplot(2,1,1)
semilogx(omega,G_prime)
xlabel('Frequency (s^-^1)')
ylabel('G''')
title('Problem 16 Graph')

subplot(2,1,2)
semilogx(omega,G_doubleprime)
xlabel('Frequency (s^-^1)')
ylabel('G''''')




%% Problem 17
% disp('Problem 17')
f_0 = 12;   % N/kg
omega_n = 10;   % rad/s
omega = 12;     % rad/s
t = 0:0.01:10;

x = 2*f_0/(omega_n^2 - omega^2)*sin((omega_n - omega)/2*t) ...
    .*sin((omega_n - omega)/2*t);

figure()
plot(t,x)
xlabel('Time (s)')
ylabel('Position (m)')
title('Problem 17 Graph')




%% Problem 18
disp('Problem 18')
R = 0.08206; % (L atm)/(mole K)
n = 1;      % moles
T = 300;    % Kelvin
a = 1.39;   % L^2 atm/mole^2
b = 0.0391; % L/mole
V = 0.08:0.02:6;    % Liters

P = n*R*T./(V-n*b) - n^2*a./V.^2;
PVRT = P.*V/(R*T);

figure()
plot(P, PVRT)
xlabel('Pressure (Pa)')
ylabel('PV/RT')
title('Problem 18 Graph')

fprintf('\tNo, the behavior of Nitrogen does not agree with the ideal gas equation\n')
fprintf('\n\n')





%% Problem 19
% disp('Problem 19')
L = 20; % meters
E = 200*10^9;   % Pa
I = 348*10^(-6); % meters^4
w = 5*10^3;     % N/m


x1 = 0:0.1:.66*L;
x2 = .66*L:0.1:L;
x = [x1 x2];

y1 = [-w*x1/(24*L*E*I).*(L*x1.^3 - 16/9*L^2*x1.^2 + 64/81*L^4)];
y2 = [-w*L/(54*E*I)*(2*x2.^3 - 6*L*x2.^2 + 40/9*L^2*x2 - 4/9*L^3)];
y = [y1 y2];

figure()
plot(x,y)
title('Problem 19 Graph')
xlabel('Distance Along Beam (m)')
ylabel('Deflection')