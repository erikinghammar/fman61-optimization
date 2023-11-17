function beta = barrier(x, A, b)
%BARRIER Returns the value of a barrier function.
    v = A*x-b;
    beta = - sum(1./(v)) + realmax*any(v>0);
end

