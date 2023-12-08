function [N_eval, varargout] = armijo(f, x, d, N_eval, varargin)
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
        switch varargin{i}
            case 'alpha'
                alpha = varargin{i+1};
            case 'epsilon'
                epsilon = varargin{i+1};
            case 'lambda'
                lambda = varargin{i+1};
            otherwise
                error(['invalid input:', variable])
        end
    end
    clear variable value;
end    

F = @(l) f(x + l*d); 

% save these to avoid unnecessary function calls.
F_0 = F(0);
N_eval = N_eval +1;
F_prim_0 = num_gradient(F,0);

% Algorithm 3, page 53 in Diehl, S. 'Optimization - A basic course'
N_eval = N_eval +1;  % this one refers to the first check below
while F(alpha*lambda) < F_0 + epsilon*F_prim_0*alpha*lambda
    % step forwards
    lambda = alpha * lambda;

    N_eval = N_eval +1;  % this one refes to the next iteration
end
N_eval = N_eval +1;  % this one refers to the first check below
while F(lambda) > F_0 + epsilon*F_prim_0*lambda
    % backtrack
    lambda = lambda/alpha;

    N_eval = N_eval +1;  % this one refes to the next iteration
end

switch nargout
    case 2
        varargout{1} = lambda;
    case 3
        varargout{1} = lambda;
        varargout{2} = F_prim_0; % this is used in wolfe_linsearch
end