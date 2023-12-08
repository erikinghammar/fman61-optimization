function[x_opt, N_eval, N_iter] = bfgs(objective_func,x0, tol, restart, printout)
%BFGS Summary of this function goes here
%   Detailed explanation goes here
%   TODO: docstring

% setup
%MAX_ITER = length(x0); % maximum number of iterations
MAX_ITER = 500;
freq = 5;

% initialization
N_eval=0;
%D_k = eye(length(x0)); % initial value for the Hessian matrix
D_k_plus = eye(length(x0));
x_opt = x0; % current best guess for optimizer.
N_iter = 0; % number of iterations
grad_k_plus = num_gradient(objective_func, x_opt);
N_eval = N_eval +2*numel(x_opt);

if printout
    lambda_k = 0;
    print_out(1, N_iter, x_opt, objective_func(x_opt), norm(grad_k_plus), N_eval, lambda_k)
end

while norm(grad_k_plus) > tol && N_iter < MAX_ITER
    grad_k = grad_k_plus;
    D_k = D_k_plus;

    % Search direction
    d_k = - D_k * grad_k;

    % line search
    [lambda_k, N_eval] = wolfe_linsearch(objective_func, x_opt, d_k, N_eval);

    %
    x_old = x_opt;
    x_opt = x_old + lambda_k*d_k;

    grad_k_plus = num_gradient(objective_func, x_opt);
    N_eval = N_eval +2*numel(x_opt);

    % p,q
    p_k = x_opt - x_old;
    q_k = grad_k_plus - grad_k;

    N_iter = N_iter + 1; % we've iterated once again.

    if printout
        % borde inte evaluera funktionen här
        print_out(0, N_iter, x_opt, objective_func(x_opt), norm(grad_k), N_eval, lambda_k)
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
    D_k_plus = D_k +                                        ...
        1/(p_k'*q_k)*(                              ...
        (1 + (q_k'*D_k*q_k)/(p_k' * q_k))       ...
        - D_k * q_k * (p_k') - p_k*(q_k')*D_k   ...
        ); % this should be the correct formula.

    if restart && mod(N_iter, freq) == 0
        D_k_plus = eye(length(x0));
    end

    if sum(D_k_plus == Inf, 'all') > 0
        disp("Stopped due to discontinuity at minimum")
        break
    end
end

end

