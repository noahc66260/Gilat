% GilatChapter8.m
% Noah Ruderman
% This script file includes all the problems from Chapter 8 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 8: Polynomials, Curve Fitting, and Interpolation
% 3.22.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
p = [0.02 -.75 0 1.25 -2];
x = -6:.1:6;
y = polyval(p,x);
figure()
plot(x, y);

fprintf('\tSee graph.\n\n');





%% Problem 2
disp('Problem 2')
p1 = [12 21 -11 -14 18 28 -4];
p2 = [4 7 -1];
[q, r] = deconv(p1, p2)





%% Problem 3
disp('Problem 3')
p1 = [4 6 -2 -5 3];
p2 = [1 4 2];
[q, r] = deconv(p1, p2)





%% Problem 4
disp('Problem 4')
m_tank = 18;
rho_steel = 7920;
d = .40;
length = .70;
r = d/2;
V = m_tank/rho_steel;

a = [-1 r];
b = [-1 length];

poly = pi*conv(b, conv(a,a));
polyeq = [0 0 0 (V - pi*r^2*length)] + poly;
x = roots(polyeq);

fprintf('\tIf the mass of the tank is 18 kg, x must be %.4f.\n\n\n', x(3));





%% Problem 5
disp('Problem 5')
length = 40;
width = 22;


disp('a.')
length_box = [-2 40];
width_box = [-2 22];
height_box = [1 0];

polyV = conv(length_box, conv(width_box, height_box))


disp('b.')
x = 0:.1:11;
V = polyV(1).*x.^3 + polyV(2).*x.^2 + polyV(3).*x + polyV(4);
figure()
plot(x,V)

fprintf('\tSee graph.\n\n')


disp('c.')
polyV1000 = polyV;
polyV1000(4) = polyV(4) - 1000;
root = roots(polyV1000);
x1 = root(2);
x2 = root(3);

fprintf('\tIf the volume of the box is 1000 cubic inches, x may be %.4f or %.4f.\n\n\n', x1, x2);


disp('d.')
max = polyder(polyV);
x = roots(max);
x1 = x(2);
Volume = polyval(polyV,x(2));

fprintf('\tThe value of x that corresponds to the largest possible volume\n')
fprintf('is %.1f inches and the largest possible volume is %.1f inches cubed.\n\n\n', ...
    x1, Volume)





%% Problem 6
disp('Problem 6')

h = .4;
V = 1.7;
r = 0.25;
polyR = [2*pi 2*pi*h/3 (r*2*pi*h/3) (-V + 2*pi*h/3*r^2)];
R = roots(polyR);
R1 = R(3);

fprintf('\tThe value of R is %.4f.\n\n\n', R1);





%% Problem 7
disp('Problem 7')
f1 = [1 -7 11 -4 -5 -2];
f2 = [0 0 0 9 -10 6];
f1plusf2 = polyadd(f1,f2,'add')
f1minusf2 = polyadd(f1,f2,'sub')





%% Problem 8
disp('Problem 8')


disp('a.')
[x,y,w] = maxormin(3, -7, 14)

disp('b.')
[x,y,w] = maxormin(-5, -11, 15)





%% Problem 9
disp('Problem 9')


disp('a.')
R = 10;
H = 30;
V_cone = 1/3*pi*R^2*H;

poly_V = -3*pi*conv([1 0], conv([1 0], [1 -10]))


disp('b.')
r = 0:.1:10;
V_cylinder = -3*pi*r.^2.*(r-10);
figure()
plot(r,V_cylinder)
xlabel('Radius of Cylinder')
ylabel('Volume of Cylinder')
title('Volume of the Cylinder vs. Radius')

fprintf('\tSee graph.\n\n');


disp('c.')
% if V = 800 in^3
V = 800;
polynom = polyadd(poly_V, V, 'sub');
r = roots(polynom);
r1 = r(1);
r2 = r(2);

fprintf('\tIf the volume of the cylinder is 800 cubic inches, the radius\n');
fprintf('may be %.4f or %.4f.\n\n', r1, r2);


disp('d.')
poly_V_der = polyder(poly_V);
max = roots(poly_V_der);
r1 = max(2);

Max_volume = polyval(poly_V, max(2));

fprintf('\tThe radius that corresponds to the maximum volume is %.4f\n', r1);
fprintf('and the maximum volume is %.4f.\n\n\n', Max_volume);





%% Problem 10
disp('Problem 10')
p = 30; % atm
T = 300; % K
a = 1.345; %L^2*atm/mole^2
b = 0.0322; %L/mole
n = 1.5;
V = waals(p, T, n, a, b)





%% Problem 11
disp('Problem 11')
x = [-6 -3.5 -2.5 -1 0 1.5 2.2 4 5.2 6 8];
y = [0.3 0.4 1.1 3.6 3.9 4.5 4.2 3.5 4.0 5.3 6.1];


disp('a.')
xplot = -6:.1:8;
polyy = polyfit(x,y,1);
yplot = polyval(polyy,xplot);
figure()
plot(x,y, 'o',xplot,yplot)
xlabel('x value')
ylabel('y value')
title('y vs. x')

fprintf('\tSee graph.\n\n');


disp('b.')
polyy = polyfit(x,y,2);
yplot = polyval(polyy,xplot);
figure()
plot(x,y, 'o',xplot,yplot)
xlabel('x value')
ylabel('y value')
title('y vs. x')

fprintf('\tSee graph.\n\n');


disp('c.')
polyy = polyfit(x,y,4);
yplot = polyval(polyy,xplot);
figure()
plot(x,y, 'o',xplot,yplot)
xlabel('x value')
ylabel('y value')
title('y vs. x')

fprintf('\tSee graph.\n\n');


disp('d.')
polyy = polyfit(x,y,10);
yplot = polyval(polyy,xplot);
figure()
plot(x,y, 'o',xplot,yplot)
xlabel('x value')
ylabel('y value')
title('y vs. x')

fprintf('\tSee graph.\n\n\n');





%% Problem 12
disp('Problem 12')
year = [1910 1930:10:2000];
population = [249 277 316 350 431 539 689 833 1014];


disp('a.')
xplot = 1910:2000;
polyy = polyfit(year, log(population), 1);
yplot = exp(polyval(polyy,xplot));
figure()
plot(xplot,yplot)
hold on;
plot(year, population, 'o')
hold off;
xlabel('Year')
ylabel('log(Population) (where population is in millions)')
title('log(Population) vs. Time in India')

pop1975 = exp(polyval(polyy,1975));

fprintf('The population in 1975 was %.1f\n\n', pop1975);


disp('b.')
polyy = polyfit(year,population,2);
yplot = polyval(polyy,xplot);
figure()
plot(year, population, 'o',xplot,yplot)
xlabel('Year')
ylabel('Population (in millions)')
title('Population vs. Time in India')

pop1975 = polyval(polyy,1975);

fprintf('\tThe population in 1975 was %.1f\n\n', pop1975);


disp('c.');
yplotlinear = interp1(year,population,xplot,'linear');
yplotspline = interp1(year,population,xplot,'spline');
figure()
plot(year, population, 'o',xplot,yplotlinear,'r',xplot,yplotspline,'g');
xlabel('Year');
ylabel('Population (in millions)');
title('Population vs. Time in India');

pop1975linear = interp1(year,population,1975,'linear');
pop1975spline = interp1(year,population,1975,'spline');

fprintf('\tThe population estimates in 1975 according to the linear and \n');
fprintf('spline fits were %.4f and %.4f million, respectively.\n\n\n', ...
    pop1975linear, pop1975spline);





%% Problem 13
disp('Problem 13')
h = [0:3:33];
D = [1.2 .91 .66 .47 .31 .19 .12 .075 .046 .029 .018 .011];


disp('a.')
figure()
plot(h,D,'o');
xlabel('Height (in km)');
ylabel('Density (in kg/m^3)');
title('Density vs. Height');

figure()
plot(log(h),D,'o');
xlabel('Log of Height (in km)');
ylabel('Density (in kg/m^3)');
title('Density vs. Log of Height');

figure()
plot(h,log(D),'o');
xlabel('Height (in km)');
ylabel('Log of Density (in kg/m^3)');
title('Log of Density vs. Height');

figure()
plot(log(h),log(D),'o');
xlabel('Log of Height (in km)');
ylabel('Log of Density (in kg/m^3)');
title('Log of Density vs. Log of Height');

fprintf('\tThe plot of log(D) vs. h gave the best fit.\n\n');


disp('b.')

figure()
plot(h,D,'o');
xlabel('Height (in km)');
ylabel('Density (in kg/m^3)');
title('Density vs. Height');

hold on

polyx = [-0.14584 0.42543];
x = 0:1:33;
y = exp(polyx(1)*x + polyx(2));
plot(x, y);
xlabel('Height (in km)');
ylabel('Density (in kg/m^3)');
title('Density vs. Height');

hold off

fprintf('\tSee graph.\n\n\n');





%% Problem 14
disp('Problem 14')

x = [0.6 2.1 3.1 5.1 6.2 7.6];
y = [0.9 9.1 24.7 58.2 105 222];
[b, m] = expofit(x,y)

figure()
plot(x,y,'o');
xlabel('x value');
ylabel('y value');
title('y vs. x with Exponential Fit');
hold on;
xplot = 0:.1:7.6;
yplot = b*exp(m*xplot);
plot(xplot, yplot)
hold off





%% Problem 15
disp('Problem 15')
v = 20:20:160;
F_D = [10 50 109 180 300 420 565 771];
rho = 1.2; % kg/m^3

polyv = polyfit(v.^2, F_D, 1);
Air_resistance = 2/rho * polyv(1);

figure()
plot(v, F_D, 'o');
xlabel('v (in m/s)');
ylabel('Drag Force (in Newtons)');
title('Drag Force vs. Velocity');
hold on;
xplot = 20:1:160;
yplot = 1/2*rho*Air_resistance*xplot.^2;
plot(xplot,yplot);
hold off;

fprintf('\tThe product CdA or the air resistance is %.4f.\n\n\n', Air_resistance)





%% Problem 16
disp('Problem 16')
t = 1:7;
P = [8.9 3.4 2.1 2.0 1.6 1.3 1.2];

figure()
plot(t,P,'o')
hold on

polynom = polyfit(sqrt(t),1./P,1);
m = polynom(1)*5
b = polynom(2)*5
xplot = 1:.1:7;
yplot = 5./(m*sqrt(xplot)+b);

plot(xplot,yplot)
hold off





%% Problem 17
disp('Problem 17')
d = [0.005 0.009 0.016 0.025 0.040 0.062 0.085 0.110];
sigma_y = [205 150 135 97 89 80 70 67];


disp('a.')

polynom = polyfit(d.^(-1/2), sigma_y, 1);
k = polynom(1)
sigma_0 = polynom(2)

yield_strength = sigma_0 + k*(0.05)^(-1/2);
fprintf('\tThe yield strength of a material with a grain size of 0.05 mm ');
fprintf('is %.0f MPa.\n\n', yield_strength);

figure()
plot(d,sigma_y,'o');
xlabel('Average Grain Diameter (in mm)');
ylabel('Yield Strength');
title('Yield Strength vs. Average Grain Diameter');
hold on
xplot = .005:.001:.110;
yplot = sigma_0 + k*xplot.^(-1/2);
plot(xplot,yplot);


disp('b.')
xplot = .005:.001:.110;
yplot = interp1(d,sigma_y,xplot,'linear');
figure()
plot(d,sigma_y,'o',xplot,yplot);
xlabel('(Average grain size)^(^-^1^/^2^) (in mm)');
ylabel('Yield Strength');
title('Yield Strength vs. (Average Grain Diameter)^(^-^1^/^2^)');

yield_strength = interp1(d,sigma_y,0.05,'linear');

fprintf('\tThe yield strength of a material with a grain size of 0.05 mm ');
fprintf('is %.0f MPa.\n\n', yield_strength);


disp('c.');
figure()
yplot = interp1(d, sigma_y, xplot,'cubic');
plot(d,sigma_y,'o',xplot,yplot);
xlabel('Average Grain Diameter (in mm)');
ylabel('Yield Strength');
title('Yield Strength vs. Average Grain Diameter');

yield_strength = interp1(d, sigma_y,0.05,'cubic');

fprintf('\tThe yield strength of a material with a grain size of 0.05 mm ');
fprintf('is %.0f MPa.\n\n\n', yield_strength);





%% Problem 18
disp('Problem 18')
V = [.75 .65 .55 .45 .35];
T = [25 37 45 56 65];
P = [1.63 1.96 2.37 3.00 3.96];
n = 0.05;

figure()
plot((T+273)./P, V, 'o')
polynom = polyfit((T+273)./P, V, 1);
R = polynom(1)/n





%% Problem 19
disp('Problem 19')

T = [-20 0 40 100 200 300 400 500 1000];
mu = [1.63 1.71 1.87 2.17 2.53 2.98 3.32 3.64 5.04];
T = T + 273;

polynom = polyfit(T, T.^(3/2)./mu, 1);
C = 1/polynom(1)
S = C*polynom(2)

T_C = -20:1:1000;
muplot = C*(T_C+273).^(3/2)./(T_C+273+S);

figure()
plot(T_C, muplot);
xlabel('Temperature (Celsius)');
ylabel('Viscosity (N s/m^2) (x 10^-^5)');
title('Viscosity vs. Temperature');
axis([-25 1000 1.5 5.1])
T=T-273;
hold on;
plot(T, mu,'o');
hold off





