function beta = barrier(x,A,b)
% barrier only works as long as points are within the constraints
bf = @(x) -1./x;  
constraints = A*x - b;
beta = sum(bf(constraints));
end