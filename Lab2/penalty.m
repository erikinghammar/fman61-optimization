function alfa = penalty(x, A, b)
% PENALTY - calculates the penalty value for an optimization problem
% 
% Description:
%   Returns the solution to alfa(x) = sum(max(0, g(x)) where g(x) := Ax - b
%   for a vector x, matrix A 
    alfa = sum(max(0, A*x -b).^2);
end

