% GilatChapter11.m
% Noah Ruderman
% This script file includes all the problems from Chapter 11 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 11: Symbolic Maths
% 5.25.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
syms x
S1 = (x - 4)^2 - (x + 3)^2 + 16*x - 4;
S2 = x^3 - 6*x^2 - x + 30;


disp('a.')
answer_1a = expand(S1 * S2);
fprintf('\t') 
disp(answer_1a)
fprintf('\b\b')


disp('b.')
answer_1b = simplify(S1 / S2);
fprintf('\t') 
disp(answer_1b)
fprintf('\b\b')


disp('c.')
answer_1c = simplify(S1 + S2);
fprintf('\t') 
disp(answer_1c)
fprintf('\b\b')


disp('d.')
answer_1d = subs(answer_1c,x,2);
fprintf('\t');
disp(answer_1d)
fprintf('\b\b')
fprintf('\n\n\n')





%% Problem 2
disp('Problem 2')
syms y x
S1 = (sqrt(3) + x)^2 - 2*(sqrt(3)*x + x/2 + x^2/2);
S2 = x^2 + 3*x + 9;


disp('a.')
answer_2a = expand(S1 * S2);
fprintf('\t') 
disp(answer_2a)
fprintf('\b\b')


disp('b.')
answer_2b = simplify(S1 / S2);
fprintf('\t') 
disp(answer_2b)
fprintf('\b\b')


disp('c.')
answer_2c = simplify(S1 + S2);
fprintf('\t') 
disp(answer_2c)
fprintf('\b\b')


disp('d.')
answer_2d = subs(answer_2c,x,4);
disp(answer_2d)
fprintf('\b\b')
fprintf('\n\n\n')





%% Problem 3
disp('Problem 3')
syms u
Q = u^3 + 2*u^2 - 25*u - 50;
R = 3*u^3 + 4*u^2 - 75*u - 100;

ans3 = simplify(Q/R);
fprintf('\t')
disp(ans3)
fprintf('\b\b')
fprintf('\n\n')





%% Problem 4
disp('Problem 4')
syms x
fx = x^5 - x^4 - 27*x^3 + 13*x^2 + 134*x - 120;


disp('a.')
fprintf('\t')
disp(factor(fx))
fprintf('\b\b')


disp('b.')
polynom = (x - 5)*(x + 3)*(x + 2)*(x - 4);
fprintf('\t')
disp(expand(polynom))
fprintf('\b\b')
fprintf('\n\n')





%% Problem 5
disp('Problem 5')
syms x y z


disp('a.')
a = 4*cos(x)^3 - 3*cos(x);
fprintf('\ta = ')
disp(a)
fprintf('\b\b')
fprintf('\tUsing the simplify command we get a = ')
disp(simplify(a))
fprintf('\b\b')


disp('b.')
b = 1/2*(sin(x - y) + sin(x + y));
fprintf('\tb = ')
disp(b)
fprintf('\b\b')
fprintf('\tUsing the simplify command we get b = ')
disp(simplify(b))
fprintf('\b\b')


disp('c.')
c = cos(x + y + z);
fprintf('\tc = ')
disp(c)
fprintf('\b\b')
fprintf('\tUsing the expand command we get c = \n\t')
disp(expand(c))
fprintf('\b\b')
fprintf('\n\n')





%% Problem 6
disp('Problem 6')
syms t
x = 3*t/(1 + t^3);
y = 3*t^2/(1 + t^3);


disp('a.')
fprintf('\tx^3 + y^3 = ')
disp(x^3 + y^3)
fprintf('\b\b\b')
fprintf(' = ')
disp(simplify(x^3 + y^3))
fprintf('\b\b')
fprintf('\t3*x*y = ')
disp(3*x*y)
fprintf('\b\b')
fprintf('\tSo we conclude that x^3 + y^3 = 3*x*y\n')


% disp('b.')
syms x y
plotfn = x^3 + y^3 - 3*x*y;
ezplot(plotfn, [-2 2 -3 2])
fprintf('\n\n')





%% Problem 7
disp('Problem 7')
SA = 370;   % m^2
h = 8;
syms R

S = pi*R*sqrt(R^2 + 4*h^2) + pi*(2*R)*h - SA;
ans7 = solve(S,R);
fprintf('\t')
fprintf('R = %.4f m \n', double(ans7(1)))
fprintf('\n\n')





%% Problem 8
disp('Problem 8')
syms a v b T T_0


disp('a.')
hilleq = (T + a)*(v + b) - (T_0 + a)*b;
hilleq_vmax = subs(hilleq, T, 0);
vmax = solve(hilleq_vmax, v);
fprintf('\t')
fprintf('v_max = ')
disp(vmax)
fprintf('\b\b')


disp('b.')
syms v_max
b_solve = solve('v_max - (b*T_0)/a = 0', b);
hilleq = subs(hilleq,b,b_solve);
fprintf('\t')
fprintf('v = ')
disp(collect(simplify(solve(hilleq,v))))
fprintf('\b\b')
fprintf('\n\n')





%% Problem 9
disp('Problem 9')
syms x y
S1 = (x - 2)^2 + (y - 3)^2 - 16;
S2 = x^2 + y^2 - 25;

figure()
ezplot(S1)
hold on
ezplot(S2)
hold off

[x y] = solve(S1, S2);
x = double(x);
y = double(y);
fprintf('\t')
fprintf('The coordinates interesect the the points (%.4f, %.4f), and \n', ...
    x(1), y(1))
fprintf('\t')
fprintf('(%.4f, %.4f)\n', x(2), y(2))
fprintf('\n\n')





%% Problem 10
disp('Problem 10')
syms T W x F_Ax F_Ay

S1 = 12/5*T - W*x;
S2 = F_Ax - 4/5*T;
S3 = F_Ay + 3/5*T - W;


disp('a.')

[T F_Ax F_Ay] = solve(S1,S2,S3,T,F_Ax,F_Ay);
fprintf('\t')
fprintf('T = ')
disp(T)
fprintf('\b\b\t')
fprintf('F_Ax = ')
disp(F_Ax)
fprintf('\b\b\t')
fprintf('F_Ay = ')
disp(F_Ay)


disp('b.')

W = 200;
T = subs(T); F_Ax = subs(F_Ax); F_Ay = subs(F_Ay);
fprintf('\t')
fprintf('T = ')
disp(T)
fprintf('\b\b\t')
fprintf('F_Ax = ')
disp(F_Ax)
fprintf('\b\b\t')
fprintf('F_Ay = ')
disp(F_Ay)
fprintf('\n')


% disp('c.')

figure()
ezplot(T, [0 4])
hold on
ezplot(F_Ax, [0 4])
ezplot(F_Ay, [0 4])
hold off
% legend('T', 'F_A_x', 'F_A_y')





%% Problem 11
disp('Problem 11')
k = .25;
syms u


% disp('a.')
p = k*u*(1 - u)/(k + u);
ezplot(p, [0 1])


disp('b.')
dp = diff(p, 'u');
u_max = solve(dp,u);
fprintf('\t')
fprintf('p is at a maximum when u = %.4f \n', double(u_max(1)))
fprintf('\t')
fprintf('The maximum value of p is %.4f \n', double(subs(p,u,u_max(1))))
fprintf('\n\n')





%% Problem 12
% disp('Problem 12')
syms x y a b m d x0 y0

S1 = x^2/a^2 + y^2/b^2 - 1;
S2 = y -m*x - d;

y1 = solve(S1,y);
y1 = y1(1);
y2 = solve(S2,y);

Dy1 = diff(y1,'x');
m0 = subs(Dy1,x,x0);
y0 = subs(y1,x,x0);
S20 = subs(S2, {y, x, m}, {y0, x0, m0});
d0 = solve(S20, 'd');

a = 10; b = 7; x0 = 8;
S1 = subs(S1);
d0 = subs(d0);
m0 = subs(m0);
S2 = subs(S2, {m, d}, {m0, d0});

figure()
ezplot(S1, [-12 12 -12 12])
hold on
ezplot(S2, [-12 12 -12 12])
hold off





%% Problem 13
disp('Problem 13')
syms theta t
h = 5;  % km
v_x = -540/60;     % km / min
x = 100;


disp('a.')
S1 = tan(theta) - h/(x + v_x*t);
theta = solve(S1, theta);
fprintf('\t')
fprintf('theta = ')
disp(theta)
fprintf('\b\b\b radians\n')


disp('b.')
dtheta = diff(theta, 't');
fprintf('\t')
fprintf('dtheta/dt = ')
disp(simplify(dtheta))
fprintf('\b\b\b radians per minute\n')
fprintf('\n\n')


% disp('c.')
theta_d = theta*180/pi;
dtheta_d = dtheta*180/pi;

figure()
ezplot(theta_d, [0 20])
axis([0 20 -90 90])

figure()
ezplot(dtheta_d, [0 20])
axis([0 20 0 110])





%% Problem 14
disp('Problem 14')
syms x

S = sin(x)^2*cos(x)/(2 + 3*sin(x))^2;
fprintf('\t')
disp(int(S, x))
fprintf('\b\b')
fprintf('\n\n')





%% Problem 15
disp('Problem 15')
syms R y H
dV = pi*R^2*(1 - y/H)^2;

V = int(dV, y, 0, H);

fprintf('\t')
fprintf('The volume of a cone is given by ')
disp(V)
fprintf('\b\b')
fprintf('\n\n')





% %% Problem 16
% disp('Problem 16')
% syms x y a b
% 
% S1 = x^2/a^2 + y^2/b^2 - 1
% x1 = solve(S1, x)
% x1 = x1(1)
% 
% w = x1
% A = int(w, y, 0, b)
% 
% % Note: The symbolic math toolbox will not solve this problem. The integral
% % cannot be evaluated analytically by MATLAB.





%% Problem 17
disp('Problem 17')
Area = 5*15;    % inches^2
syms x k A

theta = k*x;
x = 15;
theta = subs(theta) - pi;
k = solve(theta, k);
clear x
syms x
theta = k*x;

S1 = sin(theta);
S1 = A*int(S1, x, 0, 15) - Area/2;
A = solve(S1, A);

fprintf('\t')
fprintf('The value of A is ')
disp(A)
fprintf('\b\b\b')
fprintf(' and the value of k is ')
disp(k)
fprintf('\b\b')
fprintf('\n\n')





% %% Problem 18
% disp('Problem 18')
% syms x b y H x0 ybar
% 
% y = H*sin(pi*x/b);
% Area = int(y, x, 0, b)
% 
% S1 = H*sin(pi*x0/b) - ybar
% x0 = solve(S1, x0)
% x0 = x0(2)
% 
% dt = y - ybar
% t = int(dt, x, x0, b - x0)
% S2 = t - Area/2
% ybar = solve(S2, ybar)
% 
% % Note: I don't believe that this question can be done with the symbolic
% % math toolbox.





%% Problem 19
disp('Problem 19')
syms T omega V t


disp('a.')
v = V*cos(omega*t);
T = 2*pi/omega;

v_rms = sqrt(1/T * int(v^2, t, 0, T));
fprintf('\t v_rms = ')
disp(collect(v_rms))
fprintf('\b\b')


disp('b.')
v = 2.5*cos(350*t) + 3*V;
v_rms = sqrt(1/T * int(v^2, t, 0, T));
fprintf('\t v_rms = ')
disp(simplify(v_rms))
fprintf('\b\b')
fprintf('\n\n')





% %% Problem 20
% disp('Problem 20')
% syms R x N
% 
% x = dsolve('Dx = -R*x*(N + 1 - x)', 'x(0) = N', 't')
% 
% % Note: MATLAB is unable to solve for the maximum time. This is unusual as
% % I am able to solve the equations analytically with no problems. I don't
% % think MATLAB is a good program for symbolic math.



% %% Problem 21
% disp('Problem 21')
% 
% 
% 
% 
% 
% %% Problem 22
% disp('Problem 22')
% 
% 
% 
% 
% 
% %% Problem 23
% disp('Problem 23')
% 
% 
% 
% 
% 
% %% Problem 24
% disp('Problem 24')
% 
% 
% 
% 
% 
% %% Problem 25
% disp('Problem 25')
% 
% 
% 
% 
% 
