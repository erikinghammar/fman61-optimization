function [x_opt, N_eval, N_iter, normg] = nonlinearmin(f, x0, method, tol, restart, printout)
% NONLINEARMIN minimizes a given function using DFP or BFGS quasi-Newton algorithms.
%
%   [x, N_eval, N_iter, normg] = nonlinearmin(f, x0, method, tol, restart, printout)
%
% INPUTS:
%   f         - Function handle for the objective function to be minimized.
%   x0        - Initial point chosen by the user.
%   method    - String specifying the method: 'dfp' or 'bfgs'.
%   tol       - User-defined tolerance for termination.
%   restart   - Flag indicating whether to use restart (1 for yes, 0 for no).
%   printout  - Flag indicating whether to display intermediate results (1 for yes, 0 for no).
%
% OUTPUTS:
%   x         - The point at which the minimum is estimated.
%   N_eval    - Total number of function evaluations.
%   N_iter    - Total number of iterations.
%   normg     - Norm of the gradient at the output point.
%
% EXAMPLE:
%   To minimize the function stored in the m-file func.m with starting point (1, 2, 3, 4),
%   you may use the command:
%   [x, N_eval, N_iter, normg] = nonlinearmin(@func, [1, 2, 3, 4], 'DFP', 1e-6, 1, 1);
%
% NOTE:
%   The function f should be defined as a MATLAB function or anonymous function.
x_opt=x0; N_eval=0; N_iter=0;

switch lower(method)
    case 'dfp'
        % write code to perform dfp algorithm
        [x_opt, N_eval, N_iter, grad_k_plus] = dfp(f, x0, tol, restart, printout);
    case 'bfgs' % not implemented yet (?)
        % write code to perform a bfgs-algorithm
        [x_opt, N_eval, N_iter, grad_k_plus] = bfgs(f, x0, tol, restart, printout);
    otherwise
        error("non-implemented method: %s. only 'dfp' and 'bfgs' are implemented", method)
end
normg = norm(grad_k_plus);
disp("Gradient at stopping point: " + string(normg))
disp(" ")
end
