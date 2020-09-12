clc;
close all;
clear all;
N = 141;
M = 30;
S = linspace(1,30,N);
X = 10;
T = 0.5;
r = 0.05;
sigma = 0.20;
c = 0.7;
y = log(S);
dy = y(2)-y(1);
dt = T/M;
[call, put] = blsprice(S,X,r,T,sigma);
phi = zeros(N,N);
for j = 1:N
    y(j) = (j-1)*dy;
end
for i = 1:N
    for j = 1:N
        phi(i,j) = 1/sqrt((y(i) - y(j))^2 + c^2);
    end
end
%L =  zeros(N,N);
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = -(y(i)-y(j))/(sqrt((y(i) - y(j))^2 + c^2))^3;
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = -1/(sqrt((y(i) - y(j))^2 + c^2))^3 + 3*(y(i)-y(j))^2/(sqrt((y(i) - y(j))^2 + c^2))^5;
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = max(X-S,0);
%alphao = zeros(N,1);
alphao = pinv(L)*Uo';
alpha = zeros(N,1);
Unn = zeros(N,M);
t = (1:M)*dt;
Unn(1,:) = X*exp(-r*t);
for n = 1:M
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    alpha = alphao + 0.5*(F1 + F2);
    Un = L*alpha;
    %Un = X*exp(-r*n*dt);
    for j=1:N
        Unn(j,n)=Un(j);
    end
    %alpha = L\Un;
    alphao = alpha;
end
figure(1)
plot(S,Unn(:,1),'-g',S,Unn(:,round(M/2)),'-b')
    
figure(2)
mesh(t,S,Unn)
