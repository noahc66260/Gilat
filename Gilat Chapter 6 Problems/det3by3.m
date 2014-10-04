function d3 = det3by3(A)
%det3by3.m finds the determinant of a 3 by 3 matrix.
%
% Inputs:
%   A   a matrix
% 
% Outputs:
%   d3  the determinant

B1 = zeros(2,2);
B2 = zeros(2,2);
B3 = zeros(2,2);

B1(:,:) = A(2:3,2:3);
B2(1:2,1) = A(2:3,1);
B2(1:2,2) = A(2:3,3);
B3(:,:) = A(2:3,1:2);

d3 = A(1,1)*det2by2(B1) - A(1,2)*det2by2(B2) + A(1,3)*det2by2(B3);


    function d2 = det2by2(B)
    %det2by2.m finds the determinant of a 2 by 2 matrix.
        d2 = B(1,1)*B(2,2) - B(1,2)*B(2,1);
        
    end


end

