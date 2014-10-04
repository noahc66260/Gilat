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

