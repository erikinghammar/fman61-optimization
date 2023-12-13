clear
clc

%% A selection of possible inputs
%fquad = @(x) sum((x-1).^2);
fposdef = @(x, A) x'*(A*A' + 10^-3*eye(size(A)))*x;
fconv = @(x) abs(x(1)) + abs(x(2));
fnconv = @(x) sqrt(abs(x(1))+1) + sqrt(abs(x(2))+1);

% Rosenbrock function
fros = @(x) rosenbrock(x);  % dim(x) = 2

% (pen)
fpen = @(x) exp(prod(x));  % dim(x) = 5
con1 = @(x) sum(x.^2) - 10;
con2 = @(x) x(2)*x(3) - 5*x(4)*x(5);
con3 = @(x) x(1)^3 + x(3)^3 +1;
eqcons = {con1, con2, con3};
phi = 1;
pen = imposecons(fpen, phi, "penalty", eqcons, {});
% evaluates pretty fast, 0.003 seconds

dfp = "dfp";
bfgs = "bfgs";

%% Set inputs

tol = 10^-6;
restart = 0;  % 1, 0
method = dfp;
printout = 0;

f = fnconv;
x0 = [10 10]';

% f = pen;
% x0 = [-2, 2, 2, -1, -1]';  % givet av Stefan
%x0 = [0, 0, 0.0001, -1, -1]';  % bara bfgs når rätt punkt! båda med restart

% iters = zeros(2, 9);
% dims = 2.^[0 1 2 3 4 5 6 7 8];
% for i = 1:9
%     dim = dims(i)
%     x0 = 1*ones(dim, 1);
%     A = randn(length(x0));
%     f = @(x) fposdef(x, A);
%     [~, ~, iters(1, i), ~] = nonlinearmin(f, x0, dfp, tol, restart, printout);
%     [~, ~, iters(2, i), ~] = nonlinearmin(f, x0, bfgs, tol, restart, printout);
% end
% figure
% loglog(dims, iters)
% hold on
% loglog([1 2^8], [1 2^8])
% title("Iterations until convergence for different dimensions")
% legend("DFP", "BFGS", "Expected relation")

% f = fros;
% % x0 = [0.8 0.5]';
% x0 = [1.2 0.5]';

%% Run
nonlinearmin(f, x0, method, tol, restart, printout);
%[x_opt, N_evaltot, N_itertot, lastphi] = penwrapper(f, x0, method, tol, restart, printout, eqcons, {})
