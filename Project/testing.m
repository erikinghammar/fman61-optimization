clear
clc

%% Selection of inputs
fconv = @(x) x(1)^2 + x(2)^2;
dfp = "dfp";
bfgs = "bfgs";

%% Set inputs
f = fconv;
x0 = [1 1]';
method = dfp;
tol = 10^-5;
restart = 0;
printout = 1;

%% Run
[x_opt, N_eval, N_iter, normg] = nonlinearmin(f, x0, method, tol, restart, printout);