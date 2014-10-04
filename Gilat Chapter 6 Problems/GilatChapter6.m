% GilatChapter6.m
% Noah Ruderman
% This script file includes all the problems from Chapter 6 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 6: User-Defined Functions and Function Files
% 5.23.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')


disp('a.')
fprintf('\tThe fuel efficiency of a car that travels 23 miles per gallon is \n')
fprintf('%.4f kilometers per liter.\n', mgTOkm(23))


disp('b.')
fprintf('\tThe fuel efficiency of a car that travels 50 miles per gallon is \n')
fprintf('%.4f kilometers per liter.\n', mgTOkm(50))
fprintf('\n\n')





%% Problem 2
disp('Problem 2')


disp('a.')
[cm kg] = STtoSI(71,181);
fprintf('\tThe person is %.4f cm tall and weighs %.4f kg.\n', cm, kg)
fprintf('\n\n')

% Note: The book's answers deviate because of a typo. The book's answer
% corresponds to a height of 70 inches and a weight of 175 lbs.





%% Problem 3
disp('Problem 3')


disp('a.')
fprintf('\tFor x = -2, y = %9.4f\n', function3(-2))
fprintf('\tFor x =  4, y = %9.4f\n', function3(4))
fprintf('\n\n')


% disp('b.')
x = -3:.1:5;

figure()
plot(x,function3(x))
xlabel('x'); ylabel('y');
title('Problem 3 Graph Part B')





%% Problem 4
disp('Problem 4')

fprintf('\t70 km/h is %.4f feet/s.\n', kmphTOfps(70))
fprintf('\n\n')





%% Problem 5
disp('Problem 5')


disp('a.')
fprintf('\tFor theta = pi/4, r = %.4f\n', function5(pi/4))
fprintf('\tFor theta = 5*pi/4, r = %.4f\n', function5(5*pi/4))
fprintf('\n\n')

% disp('b.')

figure()
fplot(@function5, [0 2*pi])
xlabel('theta (radians)')
ylabel('r')
title('Problem 5 Graph Part B')





%% Problem 6
disp('Problem 6')


disp('a.')
[xmax y] = maxmin(2,9,-20);
fprintf('\tThe minimum of the function  2*x^2 +  9*x - 20 is at x = %7.4f.\n', ...
    xmax)
fprintf('\tAt the minimum, the value of the function is %.4f.\n', y)


disp('b.')
[xmax y] = maxmin(-3,15,50);
fprintf('\tThe maximum of the function -3*x^2 + 15*x + 50 is at x = %7.4f.\n', ...
    xmax)
fprintf('\tAt the maximum, the value of the function is %.4f.\n', y)
fprintf('\n\n')





%% Problem 7
disp('Problem 7')
fprintf('\tThe monthly payment of a 15 year mortgage of $260,000 with\n')
fprintf('annual interest rate of 6.75 percent is $%.2f.\n', ...
    amort(260000,6.75,15))
fprintf('\n\n')





%% Problem 8
disp('Problem 8')

W = @ (r,d,gamma) gamma/4*pi^2*(2*r+d)*d^2;

% W is the weight of a ring in the shape of a torus with an inner radius 'r'
% with a diameter 'd.' gamma is the specific weight of the ring material
r = 0.6;    % in
d = 0.092;  % in
gamma = 0.696;  % lb/in^3

fprintf('\tThe weight of the ring is %.4f lbs.\n', W(r,d,gamma))
fprintf('\n\n')





%% Problem 9
disp('Problem 9')


disp('a.')
a = 10; b = 15; c = 7;
fprintf('\tFor a triangle of side lengths %d, %d, and %d, the area is %.4f\n', ...
    a, b, c, triangle(a,b,c))


disp('b.')
a = 6; b = 8; c = 10;
fprintf('\tFor a triangle of side lengths %d, %d, and %d, the area is %.4f\n', ...
    a, b, c, triangle(a,b,c))


disp('c.')
a = 200; b = 75; c = 250;
fprintf('\tFor a triangle of side lengths %d, %d, and %d, the area is %.4f\n', ...
    a, b, c, triangle(a,b,c))
fprintf('\n\n')





%% Problem 10
disp('Problem 10')


disp('a.')
A = [1.5, 2.1, 4];
B = [11, 15, 9];
n = unitvec(A,B);
fprintf('\tThe unit vector in the direction of the line that connects the \n')
fprintf('points (1.5, 2.1, 4) and (11, 15, 9) is (%.4f, %.4f, %.4f).\n', ...
    n(1), n(2), n(3))


disp('b.')
A = [-11, 3, -2];
B = [-13, -4, -5];
n = unitvec(A,B);
fprintf('\tThe unit vector in the direction of the line that connects the \n')
fprintf('points (-11, 3, -2) and (-13, -4, -5) is (%.4f, %.4f, %.4f).\n', ...
    n(1), n(2), n(3))


disp('c.')
A = [1, 0, 1];
B = [0, 1, 1];
n = unitvec(A,B);
fprintf('\tThe unit vector in the direction of the line that connects the \n')
fprintf('points (1, 0, 1) and (0, 1, 1) is (%.4f, %.4f, %.4f).\n', ...
    n(1), n(2), n(3))
fprintf('\n\n')





%% Problem 11
disp('Problem 11')


disp('a.')
r1 = 5;
th1 = 23;
r2 = 12;
th2 = 40;
[r th] = AddVecPol(r1,th1,r2,th2);

fprintf('\tr = %.4f\n', r)
fprintf('\ttheta = %.4f degrees\n', th)


disp('b.')
r1 = 6;
th1 = 80;
r2 = 15;
th2 = 125;
[r th] = AddVecPol(r1,th1,r2,th2);

fprintf('\tr = %.4f\n', r)
fprintf('\ttheta = %.4f degrees\n', th)
fprintf('\n\n')





%% Problem 12
disp('Problem 12')


disp('a.')
fprintf('\tA random number between 1 and 49: %d\n', randint(1,49))


disp('b.')
fprintf('\tA random number between -35 and -2: %d\n', randint(-35,-2))
fprintf('\n\n')





%% Problem 13
disp('Problem 13')


disp('a.')
A = [1 3 2; 6 5 4; 7 8 9];
fprintf('\tThe determinant of A is %.4f\n', det3by3(A))


disp('b.')
B = [-2.5 7 1; 5 -3 -2.6; 4 2 -1];
fprintf('\tThe determinant of B is %.4f\n', det3by3(B))
fprintf('\n\n')





%% Problem 14
disp('Problem 14')


disp('a.')
A = [8, 9, 6, 10, 9, 8, 76, 86, 91, 80];
fprintf('\tThe student''s grade is %.4f\n', fgrade(A))


disp('b.')
A = [8 10 6 9 10 9 91 71 81 85;
    5 5 6 1 8 6 59 72 66 59;
    6 8 10 4 5 9 55 65 75 78;
    7 7 8 8 9 8 83 82 81 84];

grades = fgrade(A);

fprintf('\tStudent A''s grade is %.4f\n', grades(1))
fprintf('\tStudent B''s grade is %.4f\n', grades(2))
fprintf('\tStudent C''s grade is %.4f\n', grades(3))
fprintf('\tStudent D''s grade is %.4f\n', grades(4))
fprintf('\n\n')





%% Problem 15
disp('Problem 15')
R = [50, 75, 300, 60, 500, 180, 200];

fprintf('\tThe equivalent resistance is %.4f ohms\n', req(R))
fprintf('\n\n')





%% Problem 16
disp('Problem 16')
w = 250;    % mm
h = 160;    % mm
t = 26;     % mm

yc = centroidU(w,h,t);
fprintf('\tThe y coordinate of the centroid is %.4f mm\n', yc)
fprintf('\n\n')

% Note: The answer in the book is incorrect. I have been unable to find why
% the answer in the book is incorrect and could be due to a mathematical
% mistake when calculating the centroid.





%% Problem 17
disp('Problem 17')


disp('a.')
Sxx = -190;     % MPa
Syy = 145;      % MPa
Sxy = 110;      % MPa
[Smax Smin] = princstress(Sxx,Syy,Sxy);
fprintf('\tThe minimum principle stress is %+.4f MPa\n', Smin)
fprintf('\tThe maximum principle stress is %+.4f MPa\n', Smax)


disp('b.')
Sxx = 14;     % ksi
Syy = -15;    % ksi
Sxy = 8;      % ksi
[Smax Smin] = princstress(Sxx,Syy,Sxy);
fprintf('\tThe minimum principle stress is %+.4f ksi\n', Smin)
fprintf('\tThe maximum principle stress is %+.4f ksi\n', Smax)
fprintf('\n\n')





%% Problem 18
disp('Problem 18')

w = 320; h = 180; t = 32;
fprintf('\tThe moment of inertia of the beam is %.4g mm^4\n', IxcBeam(w,h,t))
fprintf('\n\n')

% There are some typos that make the answer and question unclear. I have
% calculated the centroid for U beam with the dimensions listed above. This
% is different from the book's answer and is likely related to problem 16
% where the book had an incorrect answer for the centroid.





%% Problem 19
disp('Problem 19')
fprintf('\t');
R = input('Please enter the size of the resistor in ohms: ');
fprintf('\t');
C = input('Please enter the size of the capacitor in microFarads: ');
fprintf('\n\n')

C = C*10^(-6);  % convert from mictofarads to Farads
w = 10.^(linspace(-2,6,100));

RV = lowpass(R,C,w);

figure()
semilogx(w,RV)
xlabel('\omega (rad/s)')
ylabel('RV = |V_0/V_i|')
title('Problem 19 Graph')





%% Problem 20
disp('Problem 20')

log10w = linspace(-2,7); % we are plotting the logarithm of w
w = 10.^(log10w);

fprintf('\t');
R = input('Please enter the resistance for the resistor in ohms: ');
fprintf('\t');
L = input('Please enter the inductance for the coil in Henrys: ');
fprintf('\t');
C = input('Please enter the capacitance for the capacitor in Farads: ');
fprintf('\n\n')
%i need to check if capacitance is supposed to be in farads or microfarads

R = R*ones(1,length(w));
L = L*ones(1,length(w));
C = C*ones(1,length(w));

RV = bandpass(R,C,L,w);
semilogx(w,RV)
xlabel('\omega')
ylabel('RV')
title('Problem 20 Graph RV vs. \omega')





%% Problem 21
disp('Problem 21')


disp('a.')
T = 15;     % degrees C
RH = 40;    % percent
fprintf('\tThe dew point temperature is %.4f degrees Celsius.\n', ...
    dewpoint(T,RH))


disp('b.')
T = 35;     % degrees C
RH = 80;    % percent
fprintf('\tThe dew point temperature is %.4f degrees Celsius.\n', ...
    dewpoint(T,RH))
fprintf('\n\n')




