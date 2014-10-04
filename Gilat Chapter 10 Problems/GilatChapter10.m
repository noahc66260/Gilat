% GilatChapter10.m
% Noah Ruderman
% This script file includes all the problems from Chapter 10 of MATLAB: An
% Introduction with Applications by Amos Gilat.
% Chapter 10: Applications in Numerical Analysis
% 5.19.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
y = @(x) 10*(exp(-0.2*x) - exp(-1.5*x)) - 4;

% % Let's plot this function so we can get an estimate to our zero point
% x = 0:100;
% plot(x,y(x))
% % It looks like x = 0 is a good guess

fprintf('\tThe solution is %.4f.\n\n\n', fzero(y, 0))





%% Problem 2
disp('Problem 2')
y = @(x) exp(-0.2*x).*cos(2*x) - 0.15*x.^2 + 1;

% % Let's plot our function to estimate where our zeros occur on the x axis
% x = 0:.1:4;
% plot(x,y(x))
% % Our guesses will be 1.5, 2, and 3

fprintf('\tThe positive roots are at x = %.4f, x = %.4f, and x = %.4f.\n\n\n', ...
    fzero(y,1.5), fzero(y,2), fzero(y,3))





%% Problem 3
disp('Problem 3')
y = @(x) x.^3 - 7*x.^2.*cos(3*x) + 3;

% % Let's plot this function to figure out where the positive roots are
% x = 0:.01:6.5;
% plot(x,y(x))
% % Local maxima are around 1.5, 3, 5.5. Local min are at 2, 4.2, 6.3. There
% % are six zeros

fprintf('\tThe positive roots occur at the following points:\n')
fprintf('\tx = %.4f\n', fzero(y, [1.5 2]))
fprintf('\tx = %.4f\n', fzero(y, [2 3]))
fprintf('\tx = %.4f\n', fzero(y, [3 4.2]))
fprintf('\tx = %.4f\n', fzero(y, [4.2 5.5]))
fprintf('\tx = %.4f\n', fzero(y, [5.5 6.3]))
fprintf('\tx = %.4f\n\n\n', fzero(y, [6.3 6.5]))





%% Problem 4
disp('Problem 4')
h = 10;     % m
m = 18;     % kg
mu = 0.55;   % friction coefficient
g = 9.81;   % m/s^2

y = @(x) mu*m*g*sqrt(x.^2 + h^2)./(x + mu*h) - 90;

% % Let's plot the function to see where the root will be
% x = 5:15;
% plot(x,y(x))
% % We see the zero occurs at about x = 9

fprintf('\tThe pulling force is equal to 90 N at x = %.4f meters.\n\n\n',...
    fzero(y,9))
% It appears that the book's solution is incorrect. The answer is NOT
% 7.2792





%% Problem 5
disp('Problem 5')
a = 0.18;   % meters
b = 0.06;   % meters
K = 2600;   % N/m
L_0 = sqrt(a^2 + b^2);
F = 200;    % N

L = @(x) sqrt(a^2 + (b + x).^2);
y = @(x) 2*K./(L(x)).*(L(x) - L_0).*(b + x) + K/2.*x - F;

% % Let's plot y so that we get an idea of where the zero of the function is
% x = 0:.01:.2;
% plot(x,y(x))
% % The zero is approximately x = 0.08 meters

fprintf('\tA 200 N object will hang a distance x = %.4f meters when attached\n', ...
    fzero(y,0.08))
fprintf('to the scale.\n\n\n')





%% Problem 6
disp('Problem 6')
I_s = 10^(-12);     % Amps
q = 1.6*10^(-19);   % Coulombs
k = 1.38*10^(-23);   % J/Kelvin (Boltzmann's constant)
v_s = 2;            % Volts
T = 297;            % K
R = 1000;           % Ohms

y = @(v_d) I_s*(exp(q*v_d/(k*T)) - 1) - (v_s - v_d)/R;

% % Let's plot to see where the zero is
% v_d = 0:.01:.6;
% plot(v_d, y(v_d))
% % We see that the zeros occur between x = .5 and .6

fprintf('\tv_d = %.4f Volts\n\n\n', fzero(y,[.5 .6]))





%% Problem 7
disp('Problem 7')

% we add a negative sign so that the function fminbnd will find the
% minimum of this function and thus the maximum of the regular function.
y = @(x) -(10*(exp(-0.2*x) - exp(-1.5*x)) - 4);

% % Let's plot our function to see where the maximum is
% x = 0:.1:10;
% plot(x,y(x))
% % We see that the maximum occurs between 0 and 10

[xmax ymax] = fminbnd(y, 0, 10);

fprintf('\tThe maximum of our function occurs at x = %.4f.\n', xmax)
fprintf('\tThe value of the maximum of the function is %.4f.\n\n\n', ymax)





%% Problem 8
disp('Problem 8')
% The surface area is given by pi*R1^2 + pi*(R1 + R2)*sqrt((R1 - R2)^2 +
% h^2)
% The volume is given by 1/3*pi*h*(R1^2 + R1*R2 + R2^2)
% We know R2 = 1.3*R1
% The volume is 240 cm^3
% h = 240*3/(pi*(R1^2*(2.3 + 1.3^2)))
% Our surface area is pi*R1^2 + pi*(R1*(2.3))*sqrt((R1*(0.3))^2 + (240*3/(pi*(R1^2*(2.3 + 1.3^2))))^2)
% We want to find the minimum of our surface area

SA = @(R1) pi*R1.^2 ...
    + pi*(R1.*(2.3)).*sqrt((R1.*(0.3)).^2 ... 
    + (240*3./(pi*(R1.^2*(2.3 + 1.3^2)))).^2);

% % Let's plot the SA with respect to R1 to find the minimum
% R1 = .1:.1:10;
% plot(R1, SA(R1))
% % The minimum occurs at about 4 cm

R1_min = fminbnd(SA, 3,5);
h_min = 240*3/(pi*(R1_min^2*(2.3 + 1.3^2)));

fprintf('\tFor the minimum surface area, R1 = %.4f cm and h = %.4f cm.\n\n\n', ...
    R1_min, h_min)
% Note: The answer listed in the book is incorrect. The author mistakenly
% solved the problem for V = 266 cm^3, although the problem states that V =
% 240 cm^3. 





%% Problem 9
disp('Problem 9')
h = 10;     % m
m = 18;     % kg
mu = .55;   % friction coefficient
g = 9.81;   % m/s^2

y = @(x) mu*m*g*sqrt(x.^2 + h^2)./(x + mu*h);

% % Let's plot to see the minimum of the function
% x = 0:40;
% plot(x,y(x))
% % We see there is a minimum between x = 10 and 25

[xmin Fmin] = fminbnd(y, 10, 25);

fprintf('\tThe force is smallest at x = %.4f meters. The force is %.4f N.\n\n\n', ...
    xmin, Fmin)





%% Problem 10
disp('Problem 10')
% The volume of a cone is 1/3*pi*r^2*h
% We know h < 2*R = 34 cm, r < R = 17 cm
% Based on the geometric restrictions, r = sqrt(R^2 - (h - R)^2)

R = 17;     % cm

V = @(h) 1/3*pi*(R^2 - (h - R).^2).*h;

% % Let's plot V and find a maximum
% h = 0:34;
% plot(h, V(h))

V = @(h) -1/3*pi*(R^2 - (h - R).^2).*h;
[hmax Vmax] = fminbnd(V, 0, 34);
Vmax = -Vmax
rmax = sqrt(R^2 - (hmax - R)^2);

fprintf('\tHere are the dimensions of the largest cone we can fit inside the sphere:\n')
fprintf('\tradius = %.4f cm\n', rmax)
fprintf('\theight = %.4f cm\n', hmax)
fprintf('\tvolume = %.4f cm^3 \n\n\n', Vmax)





%% Problem 11
disp('Problem 11')
c = 3.0*10^8;   % m/s
h = 6.63*10^(-34);  %J s
k = 1.38*10^(-23);  %J/K
T = 1500;   % K

% R = @(lambda) 2*pi*c^2*h./lambda.^5.*(1./(exp(h*c./(lambda*k*T)) - 1));
% lambda = linspace(.2*10^(-6), 6*10^(-6), 100);
% plot(lambda, R(lambda))

R = @(lambda) -2*pi*c^2*h./lambda.^5.*(1./(exp(h*c./(lambda*k*T)) - 1));

[lambda_max R_max] = fminbnd(R,.2*10^(-6), 6*10^(-6));
R_max = -R_max;

fprintf('\tThe maximum R is %.4g which occurs at lambda = %.4g\n\n\n', ...
    R_max, lambda_max)





%% Problem 12
disp('Problem 12')
fx = @(x) (3 + exp(.5.*x))./(.3*x.^2 + 2.5*x + 1.6);

quad(fx, 1, 6)
fprintf('\n')





%% Problem 13
disp('Problem 13')
a = 60;     % m
h = 15;     % m

fx = @(x) 2*sqrt(1 + 4*h^2/a^4*x.^2);

L = quad(fx, 0, a);
fprintf('\tThe length of the main supporting cable should be %.4f m.\n\n\n', L)





%% Problem 14
disp('Problem 14')
epsilon_0 = 8.85*10^(-12);  % C^2/N m^2
z = 5*10^(-2);  % m
R = 6*10^(-2);  % m
sigma = 300*10^(-6);    % C/m^2

E = sigma*z/(4*epsilon_0)*quad(@(r) (z^2 + r.^2).^(-3/2).*(2*r), 0, R);

fprintf('\tThe electric field is %.4g N/C\n\n\n', E)





%% Problem 15
disp('Problem 15')
R = 6371*10^3;   % m
g_0 = 9.81;     % m/s^2
m = 500;    % kg    mass of the satellite
h = 800*10^3;   % m

g = @(y) R^2./(R + y).^2 * g_0;

dU = m*quad(g, 0, h);

fprintf('\tThe change in the potential energy of the satellite is %.4g J\n\n\n', dU)





%% Problem 16
disp('Problem 16')
x = 0:30:390;
y = [42 73 101 136 143 177 205 217 212 194 190 157 141 151];

area = trapz(x,y);

fprintf('\tThe calculated area is %d square miles. This is close to the \nactual area of ', area)
fprintf('57918 square miles.\n\n\n')





%% Problem 17
disp('Problem 17')
G_inf = 5;  % ksi
c = 0.05; 
tau1 = 0.05;    % s
tau2 = 500;     % s
t = 10;     % s
% t = 100;     % s
% t = 1000;     % s

G = G_inf*(1 + c*quad(@(x) exp(-t./x)./x, tau1, tau2));

fprintf('\tFor t = %d seconds, G = %.4f\n\n\n', t, G)





%% Problem 18
disp('Problem 18')
a = 5.9065*10^9;    % km
b = 5.7208*10^9;    % km
k = sqrt(a^2 - b^2)/a;

P = 4*a*quadl(@(theta) sqrt(1 - k^2*sin(theta).^2), 0, pi/2);
v_avg = P/(248*365*24);

fprintf('\tThe perimeter of Pluto is %.4g km.\n', P)
fprintf('\tPluto has an average speed of %.4g km/hr\n\n\n', v_avg)





%% Problem 19
% disp('Problem 19')
dydt = @(x,y) x + x.*y/3;

[x y] = ode45(dydt, [1 4], 2);

figure()
plot(x,y)
xlabel('x')
ylabel('y')
title('Problem 19 Graph')





%% Problem 20
% disp('Problem 20')
dydt = @(x,y) 0.6*x.*sqrt(y) + 0.5*y.*sqrt(x);

[x y] = ode45(dydt, [0 4], 2.5);

figure()
plot(x,y)
xlabel('x')
ylabel('y')
title('Problem 20 Graph')
% Note: The book's solution is incorrect because it graphs for 
% dydt = @(x,y) -0.6*x.*sqrt(y) + 0.5*y.*sqrt(x) when the question states
% that dydt = @(x,y) 0.6*x.*sqrt(y) + 0.5*y.*sqrt(x)





%% Problem 21
disp('Problem 21')
g = 9.81;   % m/s^2
r_h = 0.025; % m

dydt = @(t, y) -sqrt(2*g*y)*r_h^2./(2 - 0.5*y).^2;

[t y] = ode45(dydt, [0 2900], 2);

time1 = t(y > .1);
time2 = t(y < .1);

time = (time2(end) - time1(end)) / 2;

figure()
plot(t,y)
xlabel('Time (s)')
ylabel('y (m)')
title('Problem 21 Graph')

fprintf('\tThe time that y = .1 is approximately %.4f seconds\n\n\n', time)





%% Problem 22
% disp('Problem 22')
C = 10^4; 
N_c = 10^4; 
r = 10^4;
R1 = 0.55;  % days^(-1)
R2 = 0.58;  % days^(-1)

dN1dt1 = @(t1, N1) R1*N1.*(1 - N1/C) - r*N1.^2./(N_c^2 + N1.^2);
dN2dt2 = @(t2, N2) R2*N2.*(1 - N2/C) - r*N2.^2./(N_c^2 + N2.^2);

[t1 N1] = ode45(dN1dt1, [0 50], 10000);
[t2 N2] = ode45(dN2dt2, [0 50], 10000);

figure()
plot(t1, N1)
title('Problem 22 Graph, R = 0.55 days^-^1')
hold on
plot(t2, N2,'--')
xlabel('Time (days)')
ylabel('Insect Number (N)')
legend('R = 0.55 days^-^1', 'R = 0.58 days^-^1')
hold off


% The book's solution is incorrect because it plots N vs. t for C = 10^5.
% However, the question states that C = 10^4, which changes the graph
% significantly. No comment can be made on the outbreak model for this
% graph because with C = 10^4 there is no outbreak. 





%% Problem 23
disp('Problem 23')
a = 5;  % lb^(1/3)
b = 2;  % day^(-1)

dwdt = @(t,w) a*w.^(2/3) - b*w;

[t w] = ode45(dwdt, [0 10], .5);

figure()
plot(t,w)
xlabel('Time (days)')
ylabel('Weight (lbs)')
title('Problem 23 Graph')

fprintf('\tThe maximum weight is %.4f pounds.\n\n\n', max(w))





%% Problem 24
% disp('Problem 24')
dvdt = @(t,v) -0.0035*v.^2 - 3;

[t v] = ode45(dvdt, [0 14.1987], 300);

figure()
subplot(2,1,1)
plot(t,v)
xlabel('Time (hrs)')
ylabel('Velocity (km/hr)')
title('Problem 24 Graph Part A')

x = zeros(size(t));
for i = 1:length(t)
    x(i) = trapz(t(1:i),v(1:i));
end

subplot(2,1,2)
plot(t,x)
xlabel('Time (hrs)')
ylabel('Distance (km)')
title('Problem 24 Graph Part B')





%% Problem 25
% disp('Problem 25')
R = 50;     % ohms
C = 0.001;  % F

%disp('a.')
v_s = 12;   % V
dv_cdt = @(t,v_c) 1/(R*C)*(v_s - v_c);

[t v_c] = ode45(dv_cdt, [0 .2], 0);

figure()
subplot(2,2,1)
plot(t,v_c)
xlabel('Time (s)')
ylabel('Capacitor Potential (V)')
title('Problem 25 Graph Part A')



%disp('b.')

dv_cdt = @(t,v_c) 1/(R*C)*(12*sin(2*60*pi*t) - v_c);

[t v_c] = ode45(dv_cdt, [0 .2], 0);

subplot(2,2,2)
plot(t,v_c)
xlabel('Time (s)')
ylabel('Capacitor Potential (V)')
title('Problem 25 Graph Part B')



%disp('c.')
v_s = 12;   % V
dv_cdt = @(t,v_c) 1/(R*C)*(v_s - v_c);

[t1 v_c1] = ode45(dv_cdt, [0 .1], 0);

v_s = 0;   % V
dv_cdt = @(t,v_c) 1/(R*C)*(v_s - v_c);

[t2 v_c2] = ode45(dv_cdt, [.1 .2], v_c1(end));

t = [t1; t2];
v_c = [v_c1; v_c2];

subplot(2,2,3)
plot(t,v_c)
xlabel('Time (s)')
ylabel('Capacitor Potential (V)')
title('Problem 25 Graph Part C')




