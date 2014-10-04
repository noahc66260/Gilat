function REQ = req(R)
%REQ.m calculates the equivalent resistance when a set of resistors is
%connected in parallel.
%
% Input:
%   R   A vector of the resistances of the resistors
%
% Output:
%   REQ     The equivalent resistance

if isempty(R)
    fprintf('There are no resistors to connect in parallel\n\n\n')
    return
end

REQ = sum(1./R);

end

