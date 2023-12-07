function [x_opt, N_eval, N_iter] = dfp(objective_func,x0, tol, restart, printout)
%DFP Summary of this function goes here
%   Detailed explanation goes here
%   TODO: write docstring
%   TODO: implement restarting condition.

% varför detta?
% MAX_ITER = length(x0)+2; % maximum number of iterations
MAX_ITER = 20; % tillfälligt

% D_k = eye(length(x0)); % initial value for the Hessian matrix
% Kasper: ändrat zeros nedan till eye, annars blir det ingen uppdatering
% och algoritmen krashar
D_k_plus = eye(length(x0));
x_opt = x0; % current best guess for optimizer.
N_iter = 0; % number of iterations
grad_k_plus = num_gradient(objective_func,x_opt);

if printout
    % borde inte evaluera funktionen här
    N_eval = 0;
    lambda_k = 0;
    print_out(1, N_iter, x_opt, objective_func(x_opt), norm(grad_k_plus), N_eval, lambda_k)
end

while norm(grad_k_plus) > tol && N_iter < MAX_ITER
    grad_k = grad_k_plus;
    D_k = D_k_plus;
    
    % Search direction
    d_k = - D_k * grad_k;

    % line search
    lambda_k = wolfe_linsearch(objective_func, x_opt, d_k);

    %
    x_old = x_opt;
    x_opt = x_old + lambda_k*d_k;

    grad_k_plus = num_gradient(objective_func, x_opt);
    
    % p,q
    p_k = x_opt - x_old;
    q_k = grad_k_plus - grad_k;

    % update Hessian
    D_k_plus = D_k + (p_k*(p_k')) / (p_k' * q_k) - ... 
               (D_k*q_k*(q_k')*D_k) / (q_k' * D_k * q_k);

    N_iter = N_iter + 1; % we've iterated once again.

    if printout
        % borde inte evaluera funktionen här
        print_out(0, N_iter, x_opt, objective_func(x_opt), norm(grad_k), N_eval, lambda_k)
    end
end



% TODO: implement proper eval counting
N_eval = -1; % number of function evaluations.
end

