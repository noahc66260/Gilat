% GilatChapter1.m
% Noah Ruderman
% This script file includes all the problems from Chapter 1 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 1: Starting with MATLAB
% 5.22.11

close all
clear 
clc
format loose

%% Problem 1
disp('Problem 1')


disp('a.')
answer = (28.5 * 3^3 - sqrt(1500)) / (11^2 + 37.3);
fprintf('\t%.4f\n', answer)


disp('b.')
answer = (7/3)^2 * 4^3 * 18 - 6^7 / (9^3 - 652);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 2
disp('Problem 2')


disp('a.')
answer = 23*(-8 + sqrt(607)/3) + (40/8 + 4.7^2)^2;
fprintf('\t%.4f\n', answer)


disp('b.')
answer = 509^(1/3) - 4.5^2 + log(200)/1.5 + 75^(1/2);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 3
disp('Problem 3')


disp('a.')
answer = (24 + 4.5^3) / (exp(4.4) - log10(12560));
fprintf('\t%.4f\n', answer)


disp('b.')
answer = 2 / 0.036 * (sqrt(250) - 10.5)^2 / exp(-.2);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 4
disp('Problem 4')


disp('a.')
answer = cos(5*pi/6)*sin(7*pi/8)^2 + tan(pi/6*log(8)) ...
    / (sqrt(7) + 2);
fprintf('\t%.4f\n', answer)


disp('b.')
answer = cos(3*pi/5)^2 + tan(pi*log(6) / 5) / (8 * 7/2);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 5
disp('Problem 5')
x = 9.75;


disp('a.') 
answer = 4*x^3 - 14*x^2 - 6.32*x + 7.3;
fprintf('\t%.4f\n', answer)


disp('b.')
answer = exp(sqrt(3)) / nthroot(0.02 * 3.1^2, 3);
fprintf('\t%.4f\n', answer)


disp('c.')
answer = log10((x^2 - x^3)^2);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 6
disp('Problem 6')
x = 5.3; z = 7.8;


disp('a.')
answer = x*z / (x/z)^2 + 14*x^2 - 0.8*z^2;
fprintf('\t%.4f\n', answer)


disp('b.')
answer = x^2*z - z^2*x + (x/z)^2 - (z/x)^(1/2);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 7
disp('Problem 7')
a = -18.2; b = 6.42; c = a/b; d = 0.5*(c*b + 2*a);


disp('a.')
answer = d - (a + b)/c + (a+d)^2/sqrt(abs(a*b*c));
fprintf('\t%.4f\n', answer)


disp('b.')
answer = log((c-d)*(b-a)) + (a + b + c + d) / (a - b - c - d);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 8
disp('Problem 8')


disp('a.')

radius = 15;
SA_sphere = 4*pi*radius^2;
SA_cube = SA_sphere;
side_cube = sqrt(SA_cube/6);
answer = side_cube;
fprintf('\t%.4f cm\n', answer)


disp('b.')

Vol_sphere = 4/3*pi*radius^3;
Vol_cube = Vol_sphere;
side_cube = nthroot(Vol_sphere, 3);
answer = side_cube;
fprintf('\t%.4f cm\n', answer)
fprintf('\n\n')





%% Problem 9
disp('Problem 9')


disp('a.')

SA_sphere = 200;
radius = sqrt(SA_sphere / (4*pi));
Vol_sphere = 4/3*pi*radius^3;
answer = Vol_sphere;
fprintf('\t%.4f in^2\n', answer)


disp('b.')
answer = 4/3*pi*sqrt(SA_sphere / (4*pi))^3;
fprintf('\t%.4f in^2\n', answer)
fprintf('\n\n')





%% Problem 10
disp('Problem 10')
x = 7/20*pi;


disp('a.')
LHS = sin(3*x);
RHS = 3*sin(x) - 4*sin(x)^3;
if LHS == RHS
    fprintf('\tBoth sides are equal\n')
else
    fprintf('\tThe sides are not equal\n')
end


disp('b.')
LHS = sin(x/2);
RHS = sqrt((1-cos(x))/2);
if LHS == RHS
    fprintf('\tBoth sides are equal\n')
else
    fprintf('\tThe sides are not equal\n')
end
fprintf('\n\n')





%% Problem 11
disp('Problem 11')
x = 27;


disp('a.')
LHS = tand(3*x);
RHS = (3*tand(x) - tand(x)^3) / (1 - 3*tand(x)^2);
if abs(LHS - RHS) < .00001
    fprintf('\tBoth sides are equal\n')
else
    fprintf('\tThe sides are not equal\n')
end


disp('b.')
LHS = tand(x/2);
RHS = sind(x) / (1 + cosd(x));
if abs(LHS - RHS) < .00001
    fprintf('\tBoth sides are equal\n')
else
    fprintf('\tThe sides are not equal\n')
end
fprintf('\n\n')





%% Problem 12
disp('Problem 12')

alpha = 5*pi/9; beta = pi/7;
LHS = sin(alpha)*sin(beta);
RHS = 1/2*(cos(alpha - beta) - cos(alpha + beta));
if abs(LHS - RHS) < .00001
    fprintf('\tBoth sides are equal\n')
else
    fprintf('\tThe sides are not equal\n')
end
fprintf('\n\n')





%% Problem 13
disp('Problem 13')

upper_bound = 3*pi/4; lower_bound = pi/3;
answer = (1/2*upper_bound - 1/4*sin(2*upper_bound))...
    - (1/2*lower_bound - 1/4*sin(2*lower_bound));
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 14
disp('Problem 14')
a = 21; b = 45; c = 60;


disp('a.')
answer = acosd((a^2 + b^2 - c^2) / (2*a*b));
fprintf('\t%.4f\n', answer)


disp('b.')
gamma = answer;
beta = asind(b*sind(gamma)/c);
alpha = asind(a*sind(gamma)/c);
fprintf('\talpha = %g \n \tbeta = %g\n', alpha, beta)


disp('c.')
sum_angles = alpha + beta + gamma;
if sum_angles == 180
    fprintf('\tThe angles add to 180 degrees\n')
else
    fprintf('\tThe angles do not add to 180 degrees\n')
end
fprintf('\n\n')





%% Problem 15
disp('Problem 15')
a = 15; b = 35;


disp('a.')
fprintf('\tc = %.4f cm\n', sqrt(a^2 + b^2))


disp('b.')
fprintf('\talpha = %.4f degrees\n', acosd(b/sqrt(a^2 + b^2)))
fprintf('\n\n')





%% Problem 16
disp('Problem 16')

A = 2; B = -7; C = -10; x_0 = 3; y_0 = -4;
d = abs(A*x_0 + B*y_0 + C) / sqrt(A^2 + B^2);

fprintf('\td = %.4f\n', d)
fprintf('\n\n')





%% Problem 17
disp('Problem 17')

fprintf('\tThe number of cartons needed is %d\n', ceil(634/18))
fprintf('\n\n')





%% Problem 18
disp('Problem 18')
CD_price = 13.95; Book_price = 44.95;
format bank


disp('a.')
answer = 3*CD_price + 5*Book_price;
fprintf('\tThe cost of 3 CD''s and 5 books is $%.2f\n', answer)


disp('b.')

answer = (3*CD_price + 5*Book_price)*1.0575;
fprintf('\tWith the 5.75%% sales tax, the price is $%.2f\n', answer)


disp('c.')
answer = round(answer);
fprintf('\tTo the nearest dollar, the price is $%d\n', answer)
fprintf('\n\n')
format short





%% Problem 19
disp('Problem 19')

answer = factorial(12) / (factorial(5)*factorial(12-5));
fprintf('\tThe number of teams that can be selected is %d\n', answer)
fprintf('\n\n')





%% Problem 20
disp('Problem 20')


disp('a.')
answer = log10(281) / log10(5);
fprintf('\t%.4f\n', answer)


disp('b.')
answer = log10(1054) / log10(7);
fprintf('\t%.4f\n', answer)
fprintf('\n\n')





%% Problem 21
disp('Problem 21')

f_0 = 100;
half_life = 3.261; % days
% we can solve for k at the half-life of Gallium-67
k = log((f_0/2)/f_0) / 3.261; % days^-1;
f_7 = f_0*exp(k*7);
f_7_rounded = round(f_7*10)/10;

fprintf('\tk = %.4f days^-1.\n', k)
fprintf('\tAfter seven days, there will be %.4f milligrams of Gallium-67 left.\n', ...
    f_7)
fprintf('\tRounded to the nearest tenth of a milligram there are %.4f milligrams.\n', ...
    f_7_rounded)
fprintf('\n\n')





%% Problem 22
disp('Problem 22')


disp('a.')
fprintf('\tThe least common mulutiple of %d and %d is %d\n', 4, 14, lcm(4,14))


disp('b.')
fprintf('\tThe least common mulutiple of %d and %d is %d\n', 8, 42, lcm(8,42))
fprintf('\n\n')





%% Problem 23
disp('Problem 23')
E_0 = 10^4.4;
energy_difference = (10^(3/2*7.1)) / 10^(3/2*6.9);
fprintf('\tThe difference in energy between earthquakes of magnitude 7.1\n')
fprintf('and 6.9 is %.4g times as many Joules\n', energy_difference)
fprintf('\n\n')





%% Problem 24
disp('Problem 24')

rate = 0.085;
P_account1 = 20000;
P_account2 = 5000;
B_account1 = P_account1*(1 + rate)^18; % after 18 years
time = log(B_account1 / P_account2) / rate;
years = floor(time);
days = ceil((time - years)*365);
fprintf('\tIt would take %d years and %d days\n', years, days)
fprintf('\n\n')





%% Problem 25
disp('Problem 25')

A = 16.0137; B = 3096.52; C = -53.67;
p_315K = exp(A - B/(C + 315));
p_405K = exp(A - B/(C + 405));
fprintf('\tAt 315 and 405 Kelvin, the vapor pressures of toluene are \n')
fprintf('%d and %d mmHg, respectively.\n', round(p_315K), round(p_405K))
fprintf('\n\n')





%% Problem 26
disp('Problem 26')

p_0 = 20*10^-6;
p_90DB = 10^(90/20)*p_0;
p_65DB = 10^(65/20)*p_0;
fprintf('\tThe sound pressure of a truck is %g Pa which is %g times\n', ...
    p_90DB, p_90DB/p_65DB)
fprintf('larger than that of a conversation.\n')
fprintf('\n\n')




