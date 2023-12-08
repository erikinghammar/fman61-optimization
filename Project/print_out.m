function print_out(first, N_itr, x_itr, fx, n_grad, ls_fun_evals, lambda)
% PRINT_ITR Print the current iteration level in the form specified in the
% manual.
%
%   print_itr(N_itr, x_itr, fx, n_grad, ls_fun_evals, lambda)
%
% INPUTS:
%   N_itr       : the current iteration level.
%   x_itr       : the current x vector.
%   fx          : the value of the function at the point x.
%   n_grad      : the norm of the gradient of the function f at point x.
%   ls_fun_eval : ?
%   

if first
    disp("iteration    x             f(x)         norm(grad)   ls fun evals   lambda")
end

disp(string(N_itr) + blanks(13 - strlength(string(N_itr))) ...
    + string(x_itr(1)) + blanks(14 - strlength(string(x_itr(1)))) ...
    + string(fx) + blanks(13 - strlength(string(fx))) ...
    + string(n_grad) + blanks(13 - strlength(string(n_grad))) ...
    + string(ls_fun_evals) + blanks(15 - strlength(string(ls_fun_evals))) ...
    + string(lambda))
for i = 2:length(x_itr)
    disp(blanks(13) + string(x_itr(i))) 
end
