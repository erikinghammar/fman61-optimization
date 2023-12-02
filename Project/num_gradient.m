function nabla_f = num_gradient(f, x0, varargin)
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
if ~isempty(varargin)
    for i=1:2:numel(varargin)
        if strcmp(varargin{i},'h')
            h = varargin{i+1}
        end
    end
end

nabla_f = zeros(size(x0));
n = numel(x0);
for i=1:n
    % for each direction we need to reset x_plus/minus_h to x and only step
    % in the considered direction.
    x_plus_h = x0;
    x_minus_h = x0;
    x_plus_h(i) = x0(i) + h;
    x_minus_h(i) = x0(i) + h;

    % central difference
    nabla_f(i) = 0.5*(f(x_plus_h)-f(x_minus_h))/h; 
end
end

