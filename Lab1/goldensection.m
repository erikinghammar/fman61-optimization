function [X,N] = goldensection(F, a, b, tol)
    phi = (sqrt(5) - 1) / 2;
    idx = 1;
    L1 = b-a;
    L(1) = L1;

    left(idx) = a;
    right(idx) = b;

    while (L(end) / L1 >= tol)
        ml = (1-phi)*(right(idx)-left(idx)) + left(idx);
        mr = phi * (right(idx)-left(idx)) + left(idx);

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

        
        L(idx) = phi * L(idx-1);
    end
    N = length(L);
    X = [left', right', L'];
end