function [X,N] = dichotomous(F, a, b, tol)
    delta = 10^-10;
    idx = 1;
    L0 = b-a;
    L(1) = L0;

    left(idx) = a;
    right(idx) = b;

    while (L(end) / L0 >= tol)
        ml = (left(idx) + right(idx))/2 - delta;
        mr = ml + 2*delta;

        Fl = F(ml);
        Fr = F(mr);

        idx = idx + 1;
        if (Fl < Fr)
            right(idx) = mr;
            left(idx) = left(idx-1);
        else 
            left(idx) = ml;
            right(idx) = right(idx-1);
        end

        
        L(idx) = abs(left(idx) - right(idx));
    end
    N = length(L)*2;
    X = [left', right', L'];
end