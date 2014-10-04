function [ cm,kg ] = STtoSI( in,lb )
%STtoSI.m takes the height and weight of a person in inches and pounds and
%converts these into centimeters and kilograms
%   in = height in inches
%   lb = weight in pounds
%   cm = height in centimeters
%   kg = weight in kilograms

cm = in*2.54;
kg = lb*0.45359237;

end

