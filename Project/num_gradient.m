function gradient = num_gradient(func, x, varargin)
%NUM_DIFF Calculate the gradient of a function 
%   Returns the gradient of a function calculated using the central
%   difference formula for each direction of the vector.
%   
%   nabla_f = num_gradient(f, x0, ['h', H_VALUE])
%
%   Inputs:
%       f           the function to be differentiated.              (Rn->R)
%       x0          the point at which to evaluate the gradient     (Rn)
%       varargin    Optional, Can be used to set the default step of the
%                       central difference method. Default h = 1e-6.
%
%   Outputs:
%       nabla_f     the gradient of f at the point x.               (Rn)
h = 1e-6;
n = numel(x);
if ~isempty(varargin)
    for i=1:2:numel(varargin)
        if strcmp(varargin{i},'h')
            h = varargin{i+1}
        end
    end
end

% Initialize gradient vector
    gradient = zeros(size(x));
    
    % Calculate partial derivatives using finite differences
    for i = 1:n
        x_plus_h = x;
        x_minus_h = x;
        x_plus_h(i) = x_plus_h(i) + h;
        x_minus_h(i) = x_minus_h(i) - h;
        
        % Calculate the partial derivative with respect to the i-th variable
        partial_derivative = (func(x_plus_h) - func(x_minus_h)) / h * 0.5;
        
        % Update the gradient vector
        gradient(i) = partial_derivative;
    end
end
