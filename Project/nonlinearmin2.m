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

if strcmp(method, 'dfp')
    updateD = @(D_k, p_k, q_k) D_k + (p_k*(p_k')) / (p_k' * q_k) - ...
        (D_k*q_k*(q_k')*D_k) / (q_k' * D_k * q_k);
elseif strcmp(method, 'bfgs')
    updateD = @(D_k, p_k, q_k) D_k +                                        ...
        1/(p_k'*q_k)*(                              ...
        (1 + (q_k'*D_k*q_k)/(p_k' * q_k))       ...
        - D_k * q_k * (p_k') - p_k*(q_k')*D_k   ...
        ); % this should be the correct formula.
else
    error("non-implemented method: %s. only 'dfp' and 'bfgs' are implemented", method)
end

% setup
MAX_ITER = 500;
freq = 5;

% initialization
N_eval=0;
%D_k = eye(length(x0)); % initial value for the Hessian matrix
D_k_plus = eye(length(x0));
x_opt = x0; % current best guess for optimizer.
N_iter = 0; % number of iterations
grad_k_plus = num_gradient(f, x_opt);

if printout
    lambda_k = 0;
    print_out(1, N_iter, x_opt, f(x_opt), norm(grad_k_plus), N_eval, lambda_k)
end

while norm(grad_k_plus) > tol && N_iter < MAX_ITER
    grad_k = grad_k_plus;
    D_k = D_k_plus;

    % Search direction
    d_k = - D_k * grad_k;

    % line search
    [lambda_k, N_eval] = wolfe_linsearch(f, x_opt, d_k, N_eval);

    %
    x_old = x_opt;
    x_opt = x_old + lambda_k*d_k;

    grad_k_plus = num_gradient(f, x_opt);
    N_eval = N_eval +2*numel(x_opt);

    % p,q
    p_k = x_opt - x_old;
    q_k = grad_k_plus - grad_k;

    N_iter = N_iter + 1; % we've iterated once again.

    if printout
        % borde inte evaluera funktionen här
        print_out(0, N_iter, x_opt, f(x_opt), norm(grad_k), N_eval, lambda_k)
    end

    if p_k == 0
        disp("Stopped due to no change in x")
        break
    elseif q_k == 0
        disp("Stopped due to no change in gradient")
        break
    elseif p_k' * q_k == 0
        disp("Stopped due to change in gradient orthogonal to change in x")
        break
    end

    if N_iter == MAX_ITER
        disp("Maximum iterations reached")
    elseif norm(grad_k_plus) <= tol
        disp("Local minimum found")
    end

    % update Hessian. Equation found on p.76 in Diehl, S. 'Optimization - a
    % basic course'.
    D_k_plus = updateD(D_k, p_k, q_k);

    if restart && mod(N_iter, freq) == 0
        D_k_plus = eye(length(x0));
    end

    if sum(D_k_plus == Inf, 'all') > 0
        disp("Stopped due to discontinuity at minimum")
        break
    end
end

normg = norm(num_gradient(f,x_opt));
disp("Gradient at stopping point: " + string(normg))
disp(" ")
end
