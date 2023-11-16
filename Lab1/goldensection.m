function [X,N] = goldensection(F, a, b, tol)
    % GOLDENSECTION Perform a golden section search for the mininum of F
    % over [a,b]
    phi = (sqrt(5) - 1) / 2;
    idx = 1;

    X = [a , b, b-a];

    R = 1; % a factor that should be 1
    while (X(end,3) / X(1,3) * R > tol)
        ml = (1-phi)*(X(idx,2)-X(idx,1)) + X(idx,1);
        mr = phi * (X(idx,2)-X(idx,1)) + X(idx,1);

        Fl = F(ml);
        Fr = F(mr);

        idx = idx + 1;
        if (Fl < Fr)
            X(idx, :) = [X(idx-1,1),mr, abs(X(idx-1,1) - mr)]; 
        else 
            X(idx, :) = [ml, X(idx-1,2), abs(X(idx-1,2) - ml)]; 
        end
    end
    N = idx;
end