function p3 = polyadd( p1,p2,operation )
%polyadd(p1,p2,operation) takes two row vectors whose elements represent
%the coefficients of polynomials and adds or subtracts them based on the
%operation specified. Acceptable operations include 'add' or 'sub'.
%   Detailed explanation goes here

% I want to make sure that the user inputs acceptable operations
while ~strcmp(operation, 'add') && ~strcmp(operation, 'sub')
    fprintf('The operation you have entered is invalid.\n');
    operation = input('Please reenter the operation. Acceptable operations are ''add'' and ''sub'': ', 's');
end

% Let's make sure that polynomial vectors are row vectors
[rows cols] = size(p1);
message = sprintf('One or more vectors have invalid dimensions. Cannot execute program.');
if rows > 1
    fprintf('%s\n', message);
    return;
end

[rows cols] = size(p2);
if rows > 1
    fprintf('%s\n', message);
    return;
end

% Let's make sure the polynomial vectors are of the same length
if length(p1) ~= length(p2)
    
    while length(p1) > length(p2)
        p2 = [0 p2];
    end
    
    while length(p2) > length(p1)
        p1 = [0 p1];
    end
end

% Now let's perform our operation
if strcmp(operation, 'add')
    p3 = p1 + p2;
else
    p3 = p1 - p2;
end



end

