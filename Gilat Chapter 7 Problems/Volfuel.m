function V = Volfuel(h)
%Volfuel.m calculates the volume of water in the shape described in
%question 22 from chapter 7 of MATLAB: An Introduction with Applications by
%Amos Gilat

r = 20;     % cm
H = 15;     % cm
L = 60;     % cm

if h < 0 || h > 55
    fprintf('Invalid input')
    return
elseif h < r
    V = L*quad(@(x) 2*((r^2 - x.^2).^(1/2) + h - r),0,(r^2 - (r - h)^2)^(1/2));
elseif h < r + H
    V = L*(pi*r^2/2 + 2*r*(h-r));
else
    V = L*(pi*r^2/2 + 2*r*H + ...
        pi*r^2/2 - quad(@(x) 2*((r^2 - x.^2).^(1/2) - (h-35)),0,(r^2 - (h-35)^2)^(1/2)));
end

    
end

