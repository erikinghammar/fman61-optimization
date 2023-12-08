clear
clc

%% Selection of inputs
fsconv = @(x) x(1)^2 + x(2)^2;
fconv = @(x) abs(x(1)) + abs(x(2));
fnconv = @(x) sqrt(abs(x(1))+1) + sqrt(abs(x(2))+1);

dfp = "dfp";
bfgs = "bfgs";

%% Set inputs
f = fsconv;
x0 = [100 100]';
% method = bfgs;
tol = 10^-6;
restart = 1;
printout = 1;

%% Run
%nonlinearmin(f, x0, dfp, tol, restart, printout);
%nonlinearmin2(f, x0, dfp, tol, restart, printout);
%nonlinearmin(f, x0, bfgs, tol, restart, printout);
nonlinearmin2(f, x0, bfgs, tol, restart, printout);
