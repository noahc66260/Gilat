function n = unitvec( A,B )
%unitvec.m find the unit vector pointing in the direction that connects the
%points A and B from A to B.
%   A = location of point A in cartesian coordinates in vector form
%   B = location of point B in cartesian coordinates in vector form
%   n = unit vector in cartesian coordinates
% We find a vector connecting A to B by subtracting the coordinates of A
% for B.

n(3) = B(3) - A(3);
n(2) = B(2) - A(2);
n(1) = B(1) - A(1);

n = n./(sqrt(n(1)^2 + n(2)^2 + n(3)^2));

end

