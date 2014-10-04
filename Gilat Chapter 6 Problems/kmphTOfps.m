function ftps = kmphTOfps(kmh)
%kmphTOfps converts kilometers per hour to feet per second
%   kmh = kilometers per hour
%   ftps = feet per second
% There are 3.2808 feet per meter

ftps = kmh * 1000 * 3.2808399 / 60 / 60;


end

