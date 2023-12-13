function [x_opt, N_evaltot, N_itertot, phi] = penwrapper(f, x0, method, tol, restart, printout, eqcons, ineqcons)

% setup
MAX_ITER = 20;

% initialization
x_opt = x0;
x_old = x_opt -1;
N_evaltot=0;
N_itertot = 0; % number of iterations
phi = 1;
iter = 0;

while iter < MAX_ITER && norm(x_opt - x_old) > tol
    x_old = x_opt;
    fwcons = imposecons(f, phi, "penalty", eqcons, ineqcons);

    [x_opt, N_eval, N_iter, ~] = nonlinearmin(fwcons, x_opt, method, tol, restart, 0);
    N_itertot = N_itertot +N_iter;
    N_evaltot = N_evaltot +N_eval;

    phi = phi *2;
    iter = iter +1;
end
phi = phi /2;

if ~printout
    clc

    disp("Optimizing the function " + func2str(f) + " from the starting point [" ...
        + num2str(x0') + "]' to the gradient tolerance " + num2str(tol) + ".")
    if restart
        disp("Restarts are used.")
    else
        disp("Restarts are not used.")
    end
    disp("Printout has been disabled.")
end

if iter == MAX_ITER
    disp("There was no convergence")
else
    disp("The solution was reached with phi = " + num2str(phi) + " in " ...
        + num2str(iter) + " phi-iterations.")
end

disp("The function value at the stopping point x is " + num2str(f(x_opt)))
disp("Final point x: [" + num2str(x_opt') + "]'.")
disp(" ")
end