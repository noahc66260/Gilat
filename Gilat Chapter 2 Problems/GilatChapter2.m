% GilatChapter2.m
% Noah Ruderman
% This script file includes all the problems from Chapter 1 of MATLAB: An
% Introduction with Applications by Amos Gilat
% Chapter 2: Creating Arrays
% 5.22.11

close all
clear
clc
format loose

%% Problem 1
disp('Problem 1')
rowvec = [6, 8*3, 81, exp(2.5), sqrt(65), sin(pi/3), 23.05];
disp(rowvec)
fprintf('\n')





%% Problem 2
disp('Problem 2')
colvec = [44; 9; log(51); 2^3; 0.1; 5*tand(25)];
disp(colvec)
fprintf('\n')





%% Problem 3
disp('Problem 3')
rowvec = 0:3:42;
disp(rowvec)
fprintf('\n')





%% Problem 4
disp('Problem 4')
colvec = [18:-4:-22]';
disp(colvec)
fprintf('\n')





%% Problem 5
disp('Problem 5')
rowvec = linspace(5,61,16);
disp(rowvec)
fprintf('\n')





%% Problem 6
disp('Problem 6')
colvec = linspace(3,-36,14);
colvec = colvec';
disp(colvec)
fprintf('\n')





%% Problem 7
disp('Problem 7')
same = zeros(1,11);
same(1,:) = 4;
disp(same)
fprintf('\n')





%% Problem 8
disp('Problem 8')
Afirst = 3:4:51
Asecond = [Afirst(1:4), Afirst(11:13)]





%% Problem 9
disp('Problem 9')
B = zeros(3,8);
B(1,:) = 0:4:28;
B(2,:) = 69:-1:62;
B(3,:) = 1.4:-.3:-.7;
disp(B)
fprintf('\n')





%% Problem 10
disp('Problem 10')
msame = zeros(3,5);
msame(:,:) = 7;
disp(msame)
fprintf('\n')





%% Problem 11
disp('Problem 11')
a = [2 -1 0 6]; b = [-5 20 12 -3]; c = [10 7 -2 1];


disp('a.')
matrix = [a; b; c];
disp(matrix)


disp('b.')
matrix = matrix';
disp(matrix)
fprintf('\n')





%% Problem 12
disp('Problem 12')


disp('a.')
a = 0:2:6;
disp(a)

disp('b.')
b = [a a];
disp(b)

disp('c.')
c = [a; a];
disp(c)

disp('d.')
d = [a' a'];
disp(d)





%% Problem 13
disp('Problem 13')
v = [2 7 -3 5 0 14 -1 10 -6 8];


disp('a.')
a = v(3:6);
disp(a)


disp('b.')
b = v([2, 4:7, 10]);
disp(b)


disp('c.')
c = v([9, 3, 1, 10]);
disp(c)


disp('d.')
d = [v([1,3 5]); v([2, 4, 6]); v([3, 6, 9])];
disp(d)





%% Problem 14
disp('Problem 14')
A = [1:5; 6:10; 11:15];


disp('a.')
va = A(1,:)


disp('b.')
vb = A(:,3)'


disp('c.')
vc = [A(2,:) A(:,4)']


disp('d.')
vd = [A(:,1)' A(:,5)']





%% Problem 15
disp('Problem 15')
B = [15:-3:3; 2:2:10; 6:6:30];


disp('a.')
ua = [B(:,2); B(:,4)]


disp('b.')
ub = B(3,:)'


disp('c.')
uc = [B(:,2)', B(:, 4)', B(:, 5)']


disp('d.')
ud = [B(:, 1)', B(1,:)]





%% Problem 16
disp('Problem 16')
A = [.1:.1:.7; 14:-2:2; 1,1,1,1,0,0,0,; 3:3:21];


disp('a.')
B = A(1:3, 1:4)


disp('b.')
C = A(2:3, :)





%% Problem 17
disp('Problem 17')
M = [6:3:21; linspace(4,4,6); 2:-1:-3; -6:2:4];


disp('a.')
A = M([1,3], [2,4]);
disp(A)


disp('b.')
B = M(:, [1, 4:6]);
disp(B)


disp('c.')
C = M([2, 3], :);
disp(C)
fprintf('\n')





%% Problem 18
disp('Problem 18')


disp('a.')
a = zeros(2,3); 
b = ones(2,3);
matrix = [a b];
disp(matrix)


disp('b.')
matrix = [b a; b a];
matrix(:, 2:5) = eye(4); 
disp(matrix)


disp('c.')
matrix = ones(4,2);
matrix(2:3, 1:2) = zeros(2);
disp(matrix)
fprintf('\n')





%% Problem 19
disp('Problem 19')
A = eye(6);
disp(A)

A(1:3, 4:6) = 3; A(5:6, 1:4) = 2;
disp(A)
fprintf('\n')





%% Problem 20
disp('Problem 20')
v = linspace(1,35,35);
newvec = reshape(v, 7, 5);
newvec = newvec';
disp(newvec)
fprintf('\n')





%% Problem 21
disp('Problem 21')
A = ones(3);
A = [A, A-A; A-A, A];
disp(A)
fprintf('\n')




