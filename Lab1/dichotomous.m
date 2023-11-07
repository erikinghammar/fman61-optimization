function [X,N] = dichotomous(F, a, b, tol)
% DICHOTOMOUS Perform a dichotomous search for the minimum value of F over
% [a,b]
    delta = 10^-10;
    idx = 1;
    X = [a, b, b-a];

    while (X(end,3) / X(1,3) >= tol)
        disp(idx)
        ml = (X(idx,1) + X(idx,2))/2 - delta;
        mr = ml + 2*delta;

        Fl = F(ml);
        Fr = F(mr);

        idx = idx + 1;
        if (Fl < Fr)
            X(idx,:) = [X(idx-1,1), mr, mr - X(idx-1,1)];
        else 
            X(idx,:) = [ml, X(idx-1,2), X(idx-1,2) - ml];
        end
    end
    N = (idx-1)*2;
end