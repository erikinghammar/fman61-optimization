function fwcons = imposecons(fwocons, phi, constype, eqcons, ineqcons) 
if strcmp(constype, "penalty")
    consf = @(x) phi*x^2;
elseif strcmp(constype, "barrier")
    consf = @(x) (x<0)*phi*(-1/x) + (x>=0)*realmax;
    if ~isempty(ineqcons)
        error("Inequality constraints cannot be combined with barrier functions.")
    end
else
    error("Only the options 'penalty' and 'barrier' are available.")
end

cons = @(x) 0;
neqcons = length(eqcons);
for i = 1:neqcons
    coni = eqcons{i};
    cons = @(x) cons(x) + consf(coni(x));
end
nineqcons = length(ineqcons);
for i = 1:nineqcons
    coni = ineqcons{i};
    cons = @(x) cons(x) + consf(coni(max(0, x)));
end

fwcons = @(x) fwocons(x) + cons(x);
end