function g = fgrade( R )
%fgrade.m calculates the final grade of students based on their homework
%assignments, midterms, and final exam
%   There are six homeworks together worth 10%. Each homework is graded on
%   a scale of 0-10. There are three midterms each worth 15%. Each midterm
%   is graded on a scale from 0-100. There is one final exam worth 45%. The
%   final exam is graded on a scale from 0-100.

%   R = a matrix in which each row are the scores of a single student. The
%   first six elements are the homework grades. The next three elements are
%   the midterm grades. The last element is the final exam grade.
%   g = a column matrix with the final grades of each student corresponding
%   to their grades in the matrix R.

[r c] = size(R);
g = zeros(r,1);

for i = 1:r
    g(i) = (sum(R(i,1:6))/60 *.1 + sum(R(i,7:9))/300*.45 ...
        + R(i,10)/100*.45)*100;
end


end

