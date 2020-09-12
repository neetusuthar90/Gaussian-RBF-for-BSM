clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5; theta = 0.5;
N = 81;
M = 100;
S = linspace(1,30,N);
dt = T/M;
y = log(S);
dy = y(2)-y(1);
%dy = 0.3;
c = 1;
[call, put] = blsprice(S,X,r,T,sigma);
for j = 1:N
    y(j) = (j-1)*dy;
end
phi = zeros(N,N);
for i = 1:N
    for j = 1:N
        phi(i,j) = sqrt((y(i) - y(j))^2 + (4*dy)^2);
    end
end
L =  zeros(N,N);
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = (y(i)-y(j))/(sqrt((y(i) - y(j))^2 + (4*dy)^2));
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 1/sqrt((y(i) - y(j))^2 + 4*dy^2)-(y(i)-y(j))^2/(sqrt((y(i) - y(j))^2 + (4*dy)^2))^3;
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
[D,E]=eig(P);
%Uo = zeros(1,N);
Uo = max(X-S,0);
K = pinv(D)*pinv(L)*Uo';
alpha = zeros(N,1);
Un = zeros(N,1);
sum = zeros(N,1);
for i = 1:100
    t = T-dt;
    for j = 1:N
        alpha = sum + K(j)*exp(E(j,j)*t)*D(:,j);
        sum = alpha;
    end
    alpha = real(alpha);
    Un = phi*alpha;
end
plot(S,Un,'k',S,put,'o')
grid on





