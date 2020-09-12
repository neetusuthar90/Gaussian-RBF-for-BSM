clc;
clear all;
close all;
X = 10;
r = 0.05;
sigma = 0.20;
T = 0.5;
t = 0;
N = 81;
S = linspace(0,30,N);
d1 = (log(S/X) + (r + 0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d2 = (log(S/X) + (r - 0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
% N1 = normcdf(d1);
% N2 = normcdf(d2);
%call_option
%V = S.*N1 - (X*exp(-r*(T-t))).*N2;
N1 = normcdf(-d1);
N2 = normcdf(-d2);
V = (X*exp(-r*(T-t))).*N2 - S.*N1;
