function [lambda] = armijo(f, x, d, varargin)
%ARMIJO Performs a line search using Armijo's algorithm
%   
%   Performs an inexact line search using Armijo's algorithm along the
%   direction d. TODO: explain more in-depth
%
%   lambda = armijo(f, x, d, ...
%           [(lambda, VALUE_LAMBDA),('alpha', VALUE_ALPHA), ('epsilon', VALUE_EPSILON)])
% 
%   Inputs:
%       f(x)        A function of x for which we are to minimize along a
%                       line d.                                     (Rn->R)
%       x           The point at which we are currently standing and from
%                       which we would like to line search in a direction.
%                                                                   (Rn)
%       d           The direction along which we line search        (Rn)
%       'lambda'    Optional, an initial guess for the lambda value. 
%                       Default is 1e-2.                            (R)
%       'alpha'     Optional, specify the value of the updated lambda
%                       factor. Default is 2.                       (R)
%       'epsilon'   Optional, specify the relative step size. Default is 0.2.
%                                                                   (R)
%   Output:
%       lambda      The value of lambda that minimizes f along line d.
%                                                                   (R)

% Initialize alpha, epsilon, lambda
alpha = 2;
epsilon = 0.2;
lambda = 1e-2;

% Check for input arguments
if ~isempty(varargin) 
    for i = 1:2:numel(varargin)
        variable = varargin{i};
        value = varargin{i+1};

        switch variable
            case 'alpha'
                alpha = value;
            case 'epsilon'
                epsilon = value;
            case 'lambda'
                lambda = value;
            otherwise
                error(['invalid input:', variable])
        end
    end
    clear variable value;
end    

F = @(l) f(x + l*d); 

% save these to avoid unnecessary function calls.
F_0 = F(0);
F_prim_0 = num_gradient(F,0);

% Algorithm 3, page 53 in Diehl, S. 'Optimization - A basic course'
while F(alpha*lambda) < F_0 + epsilon*F_prim_0*alpha*lambda
    % step forwards
    lambda = alpha * lambda;
end
while F(lambda) > F_0 + epsilon*F_prim_0*lambda
    % backtrack
    lambda = lambda/alpha;
end