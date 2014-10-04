% GilatChapter7.m
% Noah Ruderman
% This script file includes all the problems from Chapter 7 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 7: Programming in MATLAB
% 5.18.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')

disp('a.')
14 > 15/3

disp('b.')
y = 8/2 < 5*3 + 1 > 9

disp('c.')
y = 8/(2 < 5)*3 + (1 > 9)

disp('d.')
2 + 4*3 ~= 60/4 - 1
fprintf('\n')





%% Problem 2
disp('Problem 2')
a = 4;
b = 7;

disp('a.')
y = a + b >= a*b

disp('b.')
y = a + (b >= a)*b

disp('c.')
b - a < a < a/b
fprintf('\n')





%% Problem 3
disp('Problem 3')
v = [4 -2 -1 5 0 1 -3 8 2];
w = [0 2 1 -1 0 -2 4 3 2];

disp('a.')
v <= w

disp('b.')
w = v
w = [0 2 1 -1 0 -2 4 3 2];

disp('c.')
v < w + v

disp('d.')
(v < w) + v
fprintf('\n')





% Problem 4
disp('Problem 4')
v = [4 -2 -1 5 0 1 -3 8 2];
w = [0 2 1 -1 0 -2 4 3 2];

y = w(w < v)
fprintf('\n')





%% Problem 5
disp('Problem 5')

disp('a.')
-3 & 0

disp('b.')
4 < -1 & 5 > 0

disp('c.')
8 - 12 | 6 + 5 & ~-2

disp('d.')
~4 & 0 + 8*~(4 | 0)
fprintf('\n')





%% Problem 6
disp('Problem 6')
TNY = [31 26 30 33 33 39 41 41 34 33 45 42 36 39 37 45 43 36 41 37 32 32 35 42 38 33 40 37 36 51 50];   % degrees F
TANC = [37 24 28 25 21 28 46 37 36 20 24 31 34 40 43 36 34 41 42 35 38 36 35 33 42 42 37 26 20 25 31];


disp('a.')
average_NY = sum(TNY)/length(TNY);
average_ANC = sum(TANC)/length(TANC);

fprintf('\tThe average temperature in New York was %.4f degrees Fahrenheit.\n', average_NY)
fprintf('\tThe average temperature in Anchorage was %.4f degrees Fahrenheit.\n\n', average_ANC)


disp('b.')
daysbelow_NY = length(TNY(TNY < average_NY));
daysbelow_ANC = length(TANC(TANC < average_ANC));

fprintf('\tThe temperature was below average for %d days in NY.\n', daysbelow_NY)
fprintf('\tThe temperature was below average for %d days in ANC.\n\n', daysbelow_ANC)


disp('c.')
januaryDates = 1:31;
daysWarmerInANC = januaryDates(TANC > TNY);
numDays = length(daysWarmerInANC);

fprintf('\tIt was warmer in Anchorage for a total of %d days.\n', numDays)
fprintf('\tIt was warmer in Anchorage on the follwing dates in January:\n')
fprintf('\t%d ', daysWarmerInANC)
fprintf('\n\n')


disp('d.')
daysSameTemp = januaryDates(TANC == TNY);
numDays = length(daysSameTemp);

fprintf('\tIt was the same temperature in both cities for a total of %d days.\n', numDays)
fprintf('\tIt was the same temperature in both cities on the following dates in January:\n')
fprintf('\t%d ', daysSameTemp)
fprintf('\n\n')


disp('e.')
daysAboveFreezing = januaryDates(TANC > 32 & TNY > 32);
numDays = length(daysAboveFreezing);

fprintf('\tIt was above freezing in both cities for a total of %d days.\n', numDays)
fprintf('\tIt was above freezing in both cities on the following dates in January:\n')
fprintf('\t%d ', daysAboveFreezing)
fprintf('\n\n\n')





%% Problem 7
% disp('Problem 7')
x = -2:0.01:5;

% disp('a.')
fx = zeros(size(x));
for i = 1:length(x)
    if x(i) <= -1
        fx(i) = 15;
    elseif x(i) <= 1
        fx(i) = -5*x(i) + 10;
    elseif x(i) <= 3
        fx(i) = -10*x(i)^2 + 35*x(i) - 20; 
    elseif x(i) <= 4
        fx(i) = -5*x(i) + 10;
    else
        fx(i) = -10;
    end
end

figure()
subplot(2,1,1)
plot(x,fx)
title('Problem 7 Part A')
axis([-2 5 -15 20])

% disp('b.')
fx = functionFx(x);

subplot(2,1,2)
plot(x,fx)
title('Problem 7 Part B')
axis([-2 5 -15 20])





%% Problem 8
disp('Problem 8')
matrix = zeros(3,5);
for i = 1:3
    for j = 1:5
        matrix(i,j) = (i - j)/(i + j);
    end
end
matrix
fprintf('\n')





%% Problem 9
disp('Problem 9')
fprintf('\tThe quadratic equation has the form a*x^2 + b*x + c = 0.\n')
fprintf('\t');
a = input('Enter the value of a: ');
fprintf('\t');
b = input('Enter the value of b: ');
fprintf('\t');
c = input('Enter the value of c: ');

D = b^2 - 4*a*c;

if D > 0
    fprintf('\tThe equation has two roots.\n')
    root1 = (-b + sqrt(D))/(2*a);
    root2 = (-b - sqrt(D))/(2*a);
    fprintf('\tThe roots are %+.4f and %+.4f.\n', root1, root2)
elseif D == 0
    fprintf('\tThe equation has one root.\n')
    root = (-b/(2*a));
    fprintf('\tThe root is %+.4f.\n', root)
else
    fprintf('\tThe equation has no real roots.\n')
end
fprintf('\n\n')





%% Problem 10
disp('Problem 10')
m = 10;
% m = 1000;
% m = 10000;

sum = 0;
for n = 1:m
    sum = sum + 12 * (-1)^(n-1) / n^2;
end
sum = sqrt(sum);
fprintf('\tFor m = %d, our sum is %.4f.\n', m, sum)
fprintf('\tThe value of pi is %.4f.\n', pi)
fprintf('\n\n')





% Problem 11
disp('Problem 11')
x = [15 -6 0 8 -2 5 4 -10 0.5 3];

sumPos = 0;     % sum of positive numbers
sumNeg = 0;     % sum of negative numbers
for i = 1:length(x)
    if x(i) > 0
        sumPos = sumPos + x(i);
    elseif x(i) < 0
        sumNeg = sumNeg + x(i);
    end
end
fprintf('\tThe sum of the positive elements is %+4.1f.\n', sumPos)
fprintf('\tThe sum of the negative elements is %+4.1f.\n', sumNeg)
fprintf('\n\n')





%% Problem 12
disp('Problem 12')
stockReturns = [1.38 1.76 1.17 0.79 1.42 0.64 1.2 1.06 0.83 1.18];

averageReturn = Geomean(stockReturns);
fprintf('\tThe average return of the IBM stock between 1997 and 2006 was %.4f.\n', averageReturn)
fprintf('\n\n')




%% Problem 13
disp('Problem 13')

disp('a.')
fprintf('\t12! = %d \n', fact(12))
fprintf('\n')


disp('b.')
fprintf('\t0! = %d \n', fact(0))
fprintf('\n')


disp('c.')
fact(-7)
fprintf('\n')


disp('d.')
fact(6.7)
fprintf('\n')
fprintf('\n')




%% Problem 14
disp('Problem 14')

disp('a.')
fprintf('\tcos(55) = %.4f \n\n', cosTaylor(55))


disp('b.')
fprintf('\tcos(190) = %.4f \n\n', cosTaylor(190))
fprintf('\n')




%% Problem 15
disp('Problem 15')

n = 1;
found = 0;
while found == 0
    if rem(n,2) == 0
        if rem(n,7) == 0 && n^3 > 40000
            found = 1;
        else
            n = n + 2;
        end
    else
        n = n + 1;
    end
end
    
fprintf('\tThe required number is: %d\n\n\n', n)





%% Problem 16
disp('Problem 16')

disp('a.')
prime(30)


disp('b.')
prime(79)


disp('c.')
prime(52.5)


disp('d.')
prime(-20)
fprintf('\n')





%% Problem 17
disp('Problem 17')

x = randi(30, 1, 14)
downsortx = downsort(x)
fprintf('\n')





%% Problem 18
disp('Problem 18')
unsorted_matrix = randi(61,4,7)-31
sorted_matrix = matrixsort(unsorted_matrix)
fprintf('\n')





%% Problem 19
disp('Problem 19')
fprintf('\t');
age = input('Please enter your age: ');
fprintf('\t');
rhr = input('Please enter your resting heart rate: ');
fprintf('\t');
inten = input('Enter your fitness level (low, medium, or high): ', 's');

while ~strcmp(inten, 'low') && ~strcmp(inten, 'medium') && ~strcmp(inten, 'high')
    fprintf('\tNot a valid input!\n')
    fprintf('\t');
    inten = input('Enter your fitness level (low, medium, or high): ', 's');
end

if strcmp(inten, 'low')
    inten = 0.55;
elseif strcmp(inten, 'medium')
    inten = 0.65;
else
    inten = 0.8;
end

thr = int32(((220 - age) - rhr)*inten + rhr);
fprintf('\tYour training heart rate is %d.\n\n\n', thr)





%% Problem 20
disp('Problem 20')
% Note: I have chosen to write polar coordinates in the form (r, theta)
% because that is the standard for math and physics as far as I have seen
% in my undergrad studies. 

[theta radius] = CartesianToPolar(15,3);
fprintf('\t(15,3) in cartesian coordinates is (%.4f, %.4f) in polar coordinates.\n\n', radius, theta)

[theta radius] = CartesianToPolar(-7,12);
fprintf('\t(-7,12) in cartesian coordinates is (%.4f, %.4f) in polar coordinates.\n\n', radius, theta)

[theta radius] = CartesianToPolar(-17,-9);
fprintf('\t(-17,-9) in cartesian coordinates is (%.4f, %.4f) in polar coordinates.\n\n', radius, theta)

[theta radius] = CartesianToPolar(10,-6.5);
fprintf('\t(10,-6.5) in cartesian coordinates is (%.4f, %.4f) in polar coordinates.\n\n', radius, theta)
fprintf('\n')




%% Problem 21
disp('Problem 21')

fprintf('\t');
rental_period = input('Enter the rental period in days: ');

% Alert user if input is invalid and ask for new values
while rental_period < 0 || int32(rental_period) ~= rental_period
    fprintf('\tInvalid input. The rental period must be a positive integer.\n')
    fprintf('\t');
    rental_period = input('Enter the rental period in days: ');
end

fprintf('\t');
car_class = input('Enter the type of car (B, C, or D): ','s');

% Alert user if input is invalid and ask for new values
while ~strcmp(car_class, 'B') &&  ~strcmp(car_class, 'C') &&  ~strcmp(car_class, 'D')
    fprintf('\tInvalid input. The type of car must be specified by class B, C, or D.\n')
    fprintf('\t');
    car_class = input('Enter the type of car (B, C, or D): ','s');
end

if rental_period > 60
    fprintf('\tRental not available for more than 60 days.\n\n')

% Compute price for car class B
elseif strcmp(car_class,'B')
    if rental_period < 7
        price = 27*rental_period;
        fprintf('\tThe price is %d dollars. \n\n', price)
    elseif rental_period < 28
        price = 162 + 25*(rental_period - 7);
        fprintf('\tThe price is %d dollars. \n\n', price)
    else
        price = 662 + 23*(rental_period - 28);
        fprintf('\tThe price is %d dollars. \n\n', price)
    end
    
% Compute price for car class C    
elseif strcmp(car_class, 'C')
    if rental_period < 7
        price = 34*rental_period;
        fprintf('\tThe price is %d dollars. \n\n', price)
    elseif rental_period < 28
        price = 204 + 31*(rental_period - 7);
        fprintf('\tThe price is %d dollars. \n\n', price)
    else
        price = 284 + 28*(rental_period - 28);
        fprintf('\tThe price is %d dollars. \n\n', price)
    end
    
% Compute price for car class D    
else
    if rental_period < 7
        fprintf('\tClass D cars cannot be rented for less than 6 days.\n\n')
    elseif rental_period < 28
        price = 276 + 43*(rental_period - 7);
        fprintf('\tThe price is %d dollars. \n\n', price)
    else
        price = 1136 + 38*(rental_period - 28);
        fprintf('\tThe price is %d dollars. \n\n', price)
    end
end
fprintf('\n')




%% Problem 22
% disp('Problem 22')

V = zeros(1,56);
for h = 0:55
    V(h+1) = Volfuel(h);
end

h = 0:55;

figure()
plot(h,V)
xlabel('Height')
ylabel('Volume')
title('Problem 22 Graph')





%% Problem 23
% disp('Problem 23')
t = linspace(0,30,100);

x = position(t);
v = velocity(t);
a = acceleration(t);

figure()
subplot(2,2,1)
plot(t,x)
title('Problem 23 Graph Position vs. time')
xlabel('Time')
ylabel('Position')

subplot(2,2,2)
plot(t,v)
title('Velocity vs. time')
xlabel('Time')
ylabel('Velocity')

subplot(2,2,3)
plot(t,a)
title('Acceleration vs. time')
xlabel('Time')
ylabel('Acceleration')





%% Problem 24
disp('Problem 24')
format bank
bill = charge(10);
fprintf('\tThe charge is %.2f dollars.\n', bill)
fprintf('\t');
money_received = input('Enter payment as the number of the dollar bill (1, 5, or 10): ');

if money_received < bill
    fprintf('\tError. Insufficient funds.\n\n')
else
    getchange(bill, money_received);
end
format short
fprintf('\n')





%% Problem 25
% disp('Problem 25')
D_G = 150;  % mg
V_d = 50;   % L
k_a = 1.6;  % h^(-1)
k_e = 0.4;  % h^(-1)

C_p = @(t) D_G/V_d*k_a/(k_a - k_e)*(exp(-k_e*t) - exp(-k_a*t));


% disp('a.')
t = linspace(0,10,100);

figure()
subplot(2,1,1)
plot(t,C_p(t))
title('Problem 25 Part A Graph')
xlabel('Time (hours)')
ylabel('Drug Concentration')


% disp('b.')
t1 = 0:0.1:24;
t2 = 4:0.1:24;
t3 = 8:0.1:24;
t4 = 12:0.1:24;
t5 = 16:0.1:24;

C_pVector = zeros(1,length(t1));
for i = t1
    
    n = int32(i*10 + 1);
    
    if i < 4
        C_pVector(n) = C_p(i);
    elseif i < 8
        C_pVector(n) = C_p(i) + C_p(i-4);
    elseif i < 12
        C_pVector(n) = C_p(i) + C_p(i-4) + C_p(i-8);
    elseif i < 16
        C_pVector(n) = C_p(i) + C_p(i-4) + C_p(i-8) + + C_p(i-12);
    elseif i <= 24
        C_pVector(n) = C_p(i) + C_p(i-4) + C_p(i-8) + + C_p(i-12) + C_p(i-16);
    end
end

subplot(2,1,2)
plot(t1, C_pVector)
title('Problem 25 Part B Graph')
xlabel('Time (hours)')
ylabel('Drug Concentration')





%% Problem 26
disp('Problem 26')


disp('a.')
fprintf('\tThe cube root of 100 is %.4f\n\n', cubic(100))


disp('b.')
fprintf('\tThe cube root of 9261 is %.4f\n\n', cubic(9261))


disp('c.')
fprintf('\tThe cube root of -70 is %.4f\n\n', cubic(-70))
fprintf('\n')



