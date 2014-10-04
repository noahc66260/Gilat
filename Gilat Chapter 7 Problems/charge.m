function y = charge(x)
%charge.m creates a random number between 0.01 and x which is divisible by
%0.01

y = randi(x*100)/100;


end

