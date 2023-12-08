function [lambda, N_eval, F_0] = wolfe_linsearch(func, x, d, N_eval, varargin)
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
%                       Default is 0.1. 0 < epsilon <= sigma < 1.
%       'sigma',s   Optional, used to specify sigma, the tolerance. Default
%                       is 0.2.
% 0 < epsilon <= sigma < 1.
%       'alpha',a   Optional, specify alpha, the update factor. alpha > 1.
%                       Default is 2.

% TODO: Update default values once the method works.
epsilon = 0.1;
sigma = 0.2;
alpha = 5;

% Unpack the optional values.
if ~isempty(varargin) 
    for i = 1:2:numel(varargin)
        switch varargin{i}
            case 'alpha'
                if (varargin{i+1} <= 1)
                    error('alpha must be larger than 1');
                end
                alpha = varargin{i+1};
            case 'epsilon'
                epsilon = varargin{i+1};
            case 'lambda'
                lambda = varargin{i+1};
            case 'sigma'
                sigma = varargin{i+1};
        end
    end
    clear variable value;
end    

% Main method begins here.
% Method taken from page 54 of Diehl, S. 'Optimization - a basic course'

% Compute a new lambda from Armijos method.
F = @(l) func(x + l*d);
[N_eval, F_0, lambda, F_prim_0] = armijo(func,x,d,N_eval,varargin{:});

F_prim_lambda = num_gradient(F,lambda);
N_eval = N_eval +2;
a = 0;
if abs(F_prim_lambda) > -sigma * F_prim_0
    while F_prim_lambda < 0
        a = lambda;
        lambda = alpha * lambda;
        F_prim_lambda = num_gradient(F,lambda);
        N_eval = N_eval +2;
        if abs(F_prim_lambda) <= -sigma * F_prim_0
            break;
        end
    end
    b = lambda;
    lambda = (a + b)/2;
    F_prim_lambda = num_gradient(F,lambda);
    N_eval = N_eval +2;
    while abs(F_prim_lambda) > - sigma * F_prim_0
        if F_prim_lambda < 0
            a  = lambda;
        else
            b = lambda;
        end
        lambda = (a+b)/2;
        F_prim_lambda = num_gradient(F,lambda);
        N_eval = N_eval +2;
    end
end
end

