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

disp("Optimizing the function " + func2str(f) + " from the starting point [" ...
    + num2str(x0') + "]' to the gradient tolerance " + num2str(tol) + ".")
if restart
    disp("Restarts are used.")
else
    disp("Restarts are not used.")
end
if ~printout
    disp("Printout has been disabled.")
end
switch lower(method)
    case 'dfp'
        % Equation found on p.73 in Diehl, S. 'Optimization - a
        % basic course'.
        updateD = @(D_k, p_k, q_k) D_k + (p_k*(p_k')) / (p_k' * q_k) - ...
            (D_k*q_k*(q_k')*D_k) / (q_k' * D_k * q_k);
        disp("The method applied is DFP.")
    case 'bfgs'
        % Equation found on p.76 in Diehl, S. 'Optimization - a
        % basic course'.
        updateD = @(D_k, p_k, q_k) D_k +                        ...
            1/(p_k'*q_k)*(                                      ...
                (1 + (q_k'*D_k*q_k)/(p_k' * q_k))*p_k*p_k'      ...
                - D_k * q_k * (p_k') - p_k*(q_k')*D_k           ...
            );
        disp("The method applied is BFGS.")
    otherwise
        error("non-implemented method: %s. only 'dfp' and 'bfgs' are implemented", method)
end
% setup
MAX_ITER = 500;
freq = length(x0);

% initialization
N_eval=0;
D_k_plus = eye(length(x0));
x_opt = x0; % current best guess for optimizer.
x_old = x0 -1;  % just initialization
N_iter = 0; % number of iterations
grad_k_plus = num_gradient(f, x_opt);
N_eval = N_eval +2*numel(x_opt);
tic

if printout
    lambda_k = 0;
    N_eval = N_eval +1;
    print_out(N_iter, x_opt, f(x_opt), norm(grad_k_plus), N_eval, lambda_k, toc)
end

while N_iter < MAX_ITER
    grad_k = grad_k_plus;
    D_k = D_k_plus;

    % Search direction
    d_k = - D_k * grad_k;    

    % line search
    [lambda_k, N_eval, fx] = wolfe_linsearch(f, x_opt, d_k, N_eval);
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
        print_out(N_iter, x_opt, fx, norm(grad_k), N_eval, lambda_k, toc)
    end

    % glitchy stops
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

    % non-glitchy stops
    if N_iter == MAX_ITER
        disp("Maximum iterations reached")
        break
    elseif norm(grad_k_plus) <= tol
        disp("Local minimum found, gradient less than tolerance")
        break
    elseif norm(x_opt-x_old) <= tol
        disp("Local minimum found, change in x less than tolerance")
        break
    end

    D_k_plus = updateD(D_k, p_k, q_k);

    if restart && mod(N_iter, freq) == 0
        D_k_plus = eye(length(x0));
    end

    % another sort of glitchy stop
    if sum(D_k_plus == Inf, 'all') > 0
        disp("Stopped due to discontinuity at minimum")
        break
    end
end

normg = norm(grad_k_plus); 
disp("Final point x: [" + num2str(x_opt') + "]'.")
disp("Gradient norm at this x: " + string(normg))
disp(" ")
end

function N_eval = evals_add(N_eval, k)
    N_eval = N_eval +k;
end
