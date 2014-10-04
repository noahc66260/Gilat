function B = matrixsort(A)
%matrixsort.m sorts the elements of any size matrix into descending order
%across rows such that the largest value in the matrix is at position (1,1)
%and the smallest value in the matrix is at position (m,n) for an m x n
%matrix

[rowsA colsA] = size(A);
v = zeros(1,rowsA*colsA);
B = zeros(size(A));

% We will first store the values of the matrix in a row vector
for a = 1:rowsA
    for b = 1:colsA
        v((a-1)*colsA + b) = A(a,b);
    end
end

% We sort the elements of the row vector in descending order using the
% subfunction downsort(x)
v = downsort(v);

% We know store the elements of the sorted vector into the matrix B
for k = 1:rowsA
    for l = 1:colsA
         B(k,l) = v((k-1)*colsA + l);
    end
end




    function y = downsort(x)
    %downsort.m rearranges the elements of a vector x in decending order

    for i = 1:length(x)
        localMax = x(i);
        localMaxIndex = i;
        for j = i:length(x)
            if x(j) > localMax
                localMax = x(j);
                localMaxIndex = j;
            end
        end
    
        temp = x(i);
        x(i) = localMax;
        x(localMaxIndex) = temp;
    end

    y = x;


    end


end

