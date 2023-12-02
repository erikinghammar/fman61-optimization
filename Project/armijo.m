function [lambda] = armijo(f, x, d, lambda0, varargin)
%ARMIJO Performs a line search using Armijo's algorithm
%   
%   Performs an inexact line search using Armijo's algorithm along the
%   direction d. TODO: explain more in-depth
%
%   lambda = armij0(f,x,d,lambda0, ['alpha', VALUE_ALPHA, 'epsilon', VALUE_EPSILON]
% 
%   Inputs:
%       f           A function of x for which we are to minimize along a
%                       line d.
%       x           The point at which we are currently standing and from
%                       which we would like to line search in a direction.
%       d           The direction along which we want to line search.
%       lambda0     An initial guess for the lambda
%       'alpha'     Optional, specify the value of the updated lambda
%                       factor. Default is 2.
%       'epsilon'   Optional, specify the relative step size. Default is 0.2.
%
%   Output:
%       lambda

% Initialize alpha, epsilon, lambda
alpha = 2;
epsilon = 0.2;
lambda = lambda0;

% define F
F = @(l) f(x + l*d);

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
            otherwise
                error(['invalid input:', variable])
        end
    end
    clear variable value;
end    

F_0 = F(0);
F_grad_0 = gradient(F, 0);

while F(alpha*lamba) < F_0 + epsilon*F_grad_0*alpha*lambda
    lambda = alpha * lambda;
end
while F(lambda) > F_0 + epsilon*F_grad_0*lambda
    lambda = lambda/alpha;
end