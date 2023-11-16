function [X,N] = newton_wrong(F, a, b, tol)
    syms x;
    df = diff(F(x));
    d2f = diff(F(x),2);

    X = [a, b, b-a];
    idx = 2;
    lambda(1) = a;
    lambda(2) = a - double(subs(df/d2f, x, a));


    while ( ...
        abs(lambda(end) - lambda(end-1)) / (b-a) > tol ...
        && lambda(end) >= a ...
        && lambda(end) <= b ...
    )
        if (double(subs(d2f,x,lambda(idx)) < 0))
            error('d^2F/dx^2 must be positive')
        end
        lambda(idx+1) = lambda(idx) - double(subs(df/d2f,x,lambda(idx)));
        step = lambda(idx+1) - lambda(idx);
        X(idx,:) = [lambda(end)-step, lambda(end)+step, 2*step];
        idx = idx + 1;
    end
    N = idx;
end