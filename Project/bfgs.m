function[x_opt, N_eval, N_iter] = bfgs(objective_func,x0, tol, restart, printout)
%BFGS Summary of this function goes here
%   Detailed explanation goes here
%   TODO: docstring

% setup
x_opt = x0;
N_eval=0;
N_iter=0;

MAX_ITER = length(x); % maximum number of iterations

D_k = eye(length(x0)); % initial value for the Hessian matrix
D_k_plus = zeros(length(x0));
x_opt = x0; % current best guess for optimizer.
N_iter = 0; % number of iterations

grad_k = num_gradient(objective_func,x_opt);
while norm(grad_k) > tol && N_iter < MAX_ITER
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

    % update Hessian. Equation found on p.76 in Diehl, S. 'Optimization - a
    % basic course'.
    D_k_plus = D_k +                                        ...
                1/(p_k'*q_k)*(                              ...
                    (1 + (q_k'*D_k*q_k)/(p_k' * q_k))       ...
                    - D_k * q_k * (p_k') - p_k*(q_k')*D_k   ...
                ); % this should be the correct formula.

   
    N_iter = N_iter + 1; % we've iterated once again.
end

% TODO: implement proper eval counting
N_eval = -1; % number of function evaluations.

end

