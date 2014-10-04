function [lincoln washington quarter dime nickel penny] = getchange(bill, payment)
%getchange.m calculates composition of change 

lincoln = 0;
washington = 0; 
quarter = 0; 
dime = 0; 
nickel = 0; 
penny = 0;

change = round((payment - bill)*100);

while change > 0
    if change >= 500
        lincoln = lincoln + 1;
        change = change - 500;
    elseif change >= 100
        washington = washington + 1;
        change = change - 100;
    elseif change >= 25
        quarter = quarter + 1;
        change = change - 25;
    elseif change >= 10
        dime = dime + 1;
        change = change - 10;
    elseif change >= 5
        nickel = nickel + 1;
        change = change - 5;
    elseif change >= 1
        penny = penny + 1;
        change = change - 1;
    end
end


fprintf('The change consists of the following: \n')
fprintf('%-3d five dollar bills.\n', lincoln)
fprintf('%-3d one dollar bills.\n', washington)
fprintf('%-3d quarters.\n', quarter)
fprintf('%-3d dimes.\n', dime)
fprintf('%-3d nickels.\n', nickel)
fprintf('%-3d pennies.\n\n', penny)



end

