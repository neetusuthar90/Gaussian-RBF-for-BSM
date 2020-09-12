clc;close all; clear all;
X = 100; r = 0.1; sigma = 0.30;
T = 1;
N = 101;
M = 100;
y = linspace(0,6,N);
S = exp(y);
dt = T/M;
dy = y(2)-y(1);
c = 8;
[call, put] = blsprice(S,X,r,T,sigma);
for j = 1:N
    y(j) = (j-1)*dy;
end
phi = zeros(N,N);
for i = 1:N
    for j = 1:N
        phi(i,j) = exp(-(c*(y(i)-y(j)))^2);
    end
end
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = -2*c^2*(y(i)-y(j))*phi(i,j);
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = (2*c^2*phi(i,j))*(2*c*c*(y(i)-y(j))^2 - 1);
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = max(X-S,0);
%alphao = zeros(N,1);
alphao = pinv(L)*Uo';
alpha = zeros(N,1);
Unn = zeros(N,M);
t = (1:M)*dt;
%Unn(1,:) = X*exp(-r*t);
for n = 1:M
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    F3 = -dt*P*(alphao + 0.5*F2);
    F4 = -dt*P*(alphao + F3);
    alpha = alphao + (F1 + 2*F2 + 2*F3 + F4)/6;
    Un = L*alpha;
    Unn(1,n) = X*exp(-r*n*dt);
    for j=1:N
        Unn(j,n)=Un(j);
    end
    alpha = pinv(L)*Un;
    alphao = alpha;
end
% figure(1)
% plot(S,Uo)
% hold on
% plot(S,Unn(:,1),':r','linewidth',1)
% hold on
% plot(S,Unn(:,round(M/2)),'--r','linewidth',1)
% hold on
% plot(S,Unn(:,M),'-r','linewidth',1)
% hold off
% xlabel('Stock price S')
% ylabel('Option value')
% 
% legend('Option value at t = 0','Option value at t = 0.25','Option value at t = 0.5')
figure(2)
surf(t,S,Unn)
colormap cool;