clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5;
N = 81;
M = 30;
S = linspace(1,30,N);
dt = T/M;
y = log(S);
dy = y(2)-y(1);
% dy = 0.024;
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
        Ly(i,j) = (y(i)-y(j))/((abs(y(i) - y(j))^2 + (4*dy)^2));
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 1/sqrt((y(i) - y(j))^2 + (4*dy)^2)-(y(i)-y(j))^2/(sqrt((y(i) - y(j))^2 + (4*dy)^2))^3;
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
    F3 = -dt*P*(alphao + 0.5*F2);
    F4 = -dt*P*(alphao + F3);
    alpha = alphao + (F1 + 2*F2 + 2*F3 + F4)/6;
    Un = L*alpha;
    %Un(1) = X*exp(-r*n*dt);
    %alpha = L\Un;
    for j=1:N
        Unn(j,n)=Un(j);
    end
    alphao = alpha;
end
figure(1)
plot(S,Unn(:,1),'-g',S,Unn(:,round(M/2)),'-r')
pause(0.01)    
figure(2)
mesh(t,S,Unn)












