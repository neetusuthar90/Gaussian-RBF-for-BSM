clc;
clear all;
close all;
E = 10;
T = 0.5;
sigma = 0.20;
r = 0.05;
theta = 0.5;
N = 121;
Smax = 30;
M = 100;
S = linspace(1,Smax,N);
y = linspace(log(1),log(Smax),N);
dy = y(2)-y(1);
dt = T/M;
c = 4*dy;
VT = max(E - S,0);
phi = zeros(N,N);
for j = 1:N
    for i = 1:N
        phi(i,j) = sqrt(c^2 + abs(S(i) - S(j)));
    end
end
L = phi;
Ly = zeros(N,N);
Lyy = zeros(N,N);
for j = 1:N
    for i = 1:N
        Ly(i,j) = (S(i) - S(j))/phi(i,j);
    end
end
for j = 1:N
    for i = 1:N
        Lyy(i,j) = 1/phi(i,j) - (S(i) - S(j))^2/(phi(i,j))^3;
    end
end
A = zeros(N,N);
A = phi  - theta*dt*(0.5*sigma^2*S.^2*Lyy + r*S*Ly - r*phi);
H1 = phi + (1 - theta)*dt*(0.5*sigma^2*S.^2*Lyy + r*S*Ly - r*phi);
lamda = GAUSS_ELIM(phi,VT);
for i = 1:M
    t = T -i*dt;
    lamdat = 
    


