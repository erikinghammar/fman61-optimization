function [X,N] = bisection(F, a, b, tol)
%BISECTION Performs a bisect search for the minimum of F over the range
%[a,b]. 
%   Detailed explanation goes here
    idx = 1;
    X = [a, b, b-a];
    syms x;
    df = diff(F(x),x);

    while (X(end,3) / X(1,3) > tol)
        % ml = (X(idx,1) + X(idx,2))/2 - delta;
        % mr = ml + 2*delta;
        m = (X(idx,1) + X(idx,2))/2;

        % Fl = F(ml);
        % Fr = F(mr);
        % F_prim = (Fr - Fl) / (2 * delta);
        F_prim = double(subs(df, x, m));

        idx = idx + 1;
        if (F_prim > 0)
            X(idx,:) = [X(idx-1,1), m, m - X(idx-1,1)];
        else 
            X(idx,:) = [m, X(idx-1,2), X(idx-1,2) - m];
        end
    end
    N = (idx);
end