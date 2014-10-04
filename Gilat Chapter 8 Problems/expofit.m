function [ b,m ] = expofit( x,y )
%expofit(x,y) takes vectors x, y and returns the values of b and m in the
%equation y = b*e^(m*x) after fitting the data in that form.

polynom = polyfit(x,log(y),1);
b = exp(polynom(2));
m = polynom(1);

end

