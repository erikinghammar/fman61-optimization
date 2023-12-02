function Df = num_diff(f,varargin)
%NUM_DIFF Differentiate a function 
%   Differentiate a function using the central difference formula:
%   f'(x) ≈ (f(x+h) - f(x-h)) / (2h). This gets better as h->0.
%   
%   Inputs:
%       f       the function to be differentiated                   (R->R)
%       h       Optional, the size of the step in the central difference.
%                   Default value is 1e-6.
%   Outputs:
%       Df      function handle for the numerical derivative of f.  (R)

h = 1e-6;
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        variable = varargin{i};
        value = varargin{i+1};

        switch variable % this can be modified to add more variables if needed.
            case 'h'
                h = value
        end
    end
    clear variable value;
end
Df = @(x) (f(x+h) - f(x-h)) / (2*h);
end

