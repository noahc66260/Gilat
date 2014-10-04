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
