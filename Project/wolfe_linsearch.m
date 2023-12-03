function [lambda] = wolfe_linsearch(func, x, d, varargin)
%WOLF_LINSEARCH Perform a line search to satisfy the Wolfe condition
%   The Wolfe conditions creates an interval of acceptable points to be
%   used when performing an inexact line search. They are:
%       F(lambda) <= F(0) + epsilon*F'(0)*lambda;
%       abs(F'(lambda))<= -sigma*F'(0);
%   Where:
%       F(l) = f(x + l*d), f:Rn->R, x, d: Rn, l: R.
%       0 < epsilon <= sigma < 1
%   
%   This algorithm uses Armijo's method to estimate an initial lambda0 and
%   then further tries to find a value that satisfies the Wolfe conditions.
%
%   lambda = wolfe_linsearch(func, x, d, ...
%       ('lambda', LAMBDA_GUESS), ('epsilon',EPSILON_VALUE'), ...
%       ('sigma', SIGMA_VALUE), ('alpha', ALPHA_VALUE)
%   )
%
%   Inputs:
%       func        objective function
%       x           point from which to line search.
%       d           direction in which to search.
%       'lambda',l  Optional, can be used to guess an initial
%                       lambda value.
%       'epsilon',e Optional, used to specify the tolerance of the search.
%                       Default is 0.2. 0 < epsilon <= sigma < 1.
%       'sigma',s   Optional, used to specify sigma, the tolerance. Default
%                       is 0.9. 0 < epsilon <= sigma < 1.
%       'alpha',a   Optional, specify alpha, the update factor. alpha > 1.
%                       Default is 2.

% TODO: Update default values once the method works.
epsilon = 0.2;
sigma = 0.9;
alpha = 2;

% Unpack the optional values.
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
            case 'sigma'
                sigma = value;
        end
    end
    clear variable value;
end    

% Main method begins here.
% Method taken from page 54 of Diehl, S. 'Optimization - a basic course'

% Compute a new lambda from Armijos method.
F = @(l) func(x + l*d);
[lambda, F_prim_0] = armijo(func,x,d,varargin{:});
F_prim_lambda = num_gradient(F,lambda);
a = 0;
if abs(F_prim_lambda) > -sigma * F_prim_0
    while F_prim_lambda < 0
        a = lambda;
        lambda = alpha * lambda;
        F_prim_lambda = num_gradient(F,lambda);
        if abs(F_prim_lambda) <= -sigma * F_prim_0
            break;
        end
    end
    b = lambda;
    lambda = (a + b)/2;
    F_prim_lambda = num_gradient(F,lambda);
    while abs(F_prim_lambda) > - sigma * F_prim_0
        if F_prim_lambda < 0
            a  = lambda;
        else
            b = lambda;
        end
        lambda = (a+b)/2;
        F_prim_lambda = num_gradient(F,lambda);
    end
end
end

