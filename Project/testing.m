clear
clc

%% Selection of inputs
fconv = @(x) x(1)^2 + x(2)^2;
fnconv = @(x) sqrt(abs(x(1))+1) + sqrt(abs(x(2))+1);

dfp = "dfp";
bfgs = "bfgs";

%% Set inputs
f = @(x) abs(x(1)) + abs(x(2));
x0 = [1.1 1.1]';
% method = bfgs;
tol = 10^-6;
restart = 0;
printout = 1;

%% Run
nonlinearmin(f, x0, dfp, tol, restart, printout);
nonlinearmin(f, x0, bfgs, tol, restart, printout);