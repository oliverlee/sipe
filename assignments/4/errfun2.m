function [e, ymod] = errfun2(par)

global t y

% Define parameters from parameter vector
a = par(1);
b = par(2);

% Define estimated y: ymod
ymod = a*cos(b*t/2) + b*sin(a*t/5);

% Define error: e
e = y - ymod;

