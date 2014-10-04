function n = randint( a,b )
%randint.m gives a random integer number between a range of two numbers
%   a = lower range
%   b = upper range
%   n = random number between a and b

range = b - a + 1;
n = randi(range) + a - 1;


end

