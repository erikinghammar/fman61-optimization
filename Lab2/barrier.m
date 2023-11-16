function beta = barrier(x, A, b)
%BARRIER Returns the value of a barrier function.
    beta = - sum(1./(A*x-b));
end

