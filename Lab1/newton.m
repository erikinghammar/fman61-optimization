function [X, N] = newton(dF, d2F, lambda0, tol)
    idx  = 1;
    MAX_ITR = 100;
    lambda(1) = lambda0;

    while ( ...
        abs(dF(lambda(end)))>tol ...
        && (idx < MAX_ITR) ...
    )
        idx = idx+1;
        lambda(idx) = lambda(idx-1) - dF(lambda(idx-1))/d2F(lambda(idx-1));
    end
    N = idx - 1;
    X = lambda(end);
    if (isnan(X))
        error('no convergence')
    end
end