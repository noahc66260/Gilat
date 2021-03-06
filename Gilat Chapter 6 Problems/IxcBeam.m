function Ixc = IxcBeam(w,h,t)
%IxcBeam.m finds the area moment of intertia of a U beam about the axis x_0
%that passes through its centroid.
%   w = width of the object
%   h = the height of the arms of the boject
%   t = the width of the arms of the object and the height of the central
%   bar

yc = centroidU(w,h,t);

centroid_y1 = t/2;
centroid_y23 = h/2;

I_x01 = 1/12*(w - 2*t)*t^3 + (w - 2*t)*t*(yc - centroid_y1)^2;
I_x02 = 1/12*(t)*h^3 + (t*h)*(yc - centroid_y23)^2;
I_x03 = I_x02;

Ixc = I_x01 + I_x02 + I_x03;


function yc = centroidU(w,h,t)
%centroidU.m calculates the centroid along the y axis of U-shaped cross
%sectional figure shown problem 16 of Gilat
%   w = width of the object
%   h = the height of the arms of the boject
%   t = the width of the arms of the object and the height of the central
%   bar

if 2*t*(h-t) > w*t
    yc = (w-2*t-2*h)/(-4);
    % 2*t*(h-yc) = 2*t*(yc-t)+t*w 
    % -2t*yc +2th = 2tyc -2t^2 +tw
    % -4tyc = tw -2t^2 -2th
    % yc = (w-2*t-2*h)/(-4)
else
    yc = (2*t*h - 2*t^2 +w*t)/(2*w);
    % 2*t*(h-t) + (t-yc)*w = yc*w
    % 2t(h-t) = w(yc - (t-yc))
    % w(2yc - t) = 2t(h-t)
    % yc*2w - wt = 2th - 2t^2
    % yc*2w = 2t*h -2t^2 + wt
    % yc = (2*t*h - 2*t^2 +w*t)/(2*w)

end
end


end

