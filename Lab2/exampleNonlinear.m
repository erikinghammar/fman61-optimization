clear; clc; close all

%% INPUT DATA
% feasible region from Problem 1
A =  [-1, -1;1, -1;1, 2;-1, 0];
b = [-2;0;6;0];

% for compatibility, function must be written with two input arguments
% as opposed to f = @(x) x(1) + x(2) ...
f = @(x,y) x.^2 + y.^2;

x0 = [0;0]; % initial point

% plot feasible region + level curves
feasibleRegion(A,b,f);

%% OPTIMIZATION
% penalty method
mu = 1e1;
alpha = @(x) penalty(x,A,b);
F = @(x) f(x(1),x(2)) + mu*alpha(x);

% unconstrained optimization
options = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton');
xopt = fminunc(F,x0,options);
fprintf('x_optimal = (%f,%f)\n',xopt);
hold on
plot(xopt(1),xopt(2),'ro'); axis square;
hold off

%%
function [p,c] = feasibleRegion(A,b,f)
% plots the feasible region for the linear programming problem
%   minimize f(x) 
% s.t. Ax <= b
%
% INPUTS
% A: matrix of constraints (2x2 or 3x3)
% b: vector of constraints (2x1 or 3x1)
% f: anonymous function of 2 or 3 inputs, e.g. f = @(x,y) ... f = @(x,y,z)
%
% OUTPUTS
% p: fill object (2D) or alphaShape (3D) with feasible region in yellow
% f: level curves for objective function (in 2D case)

if nargin(f) ~= size(A,2)
    error('No. of arguments in function f does not match number of columns in A.');
end

comb = nchoosek(1:size(A,1),size(A,2));
v = [];
for k = 1:size(comb,1)
    Ak = A(comb(k,:),:);
    bk = b(comb(k,:));
    if rank(Ak) == size(A,2)
        xk = Ak\bk;
        if all(A*xk - b(:) <= 1e-12)
            v = [v xk];
        end
    end
end

if length(v) < 3
    error('Feasible region is not bounded.')
end

figure
switch size(A,2)
    case 2
        set(gca,'View',[0 90]);
        [~,I] = unique(angle(v(1,:) + 1i*v(2,:) - mean(v(1,:) + 1i*v(2,:))));
        p = fill(gca,v(1,I),v(2,I),'y','FaceAlpha',0.5,'HitTest','off');
    case 3
        v = unique(v.','rows');
        set(app.UIAxes,'View',[-37.5, 30]);
        p = plot(alphaShape(v(:,1),v(:,2),v(:,3)),'Parent',gca,'FaceAlpha',0.5,'FaceColor','yellow');
end

if nargin > 2 && size(A,2) == 2
    x = gca().XLim;
    y = gca().YLim;
    [x,y] = meshgrid(linspace(x(1),x(2)),linspace(y(1),y(2)));
    z = f(x,y);
    s = -1.1*min(min(z(:)),0) + 1e-6;
    t = logspace(log10(min(z(:)) + s),log10(max(z(:)) + s),20) - s;
    hold on
    c = contour(x,y,z,t);
    hold off
    colorbar; colormap winter;
    grid on
end
end