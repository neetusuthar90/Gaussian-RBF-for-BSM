clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5;
N =81;
M = 30;
S = linspace(1,30,N);
dt = T/M;
y = log(S);
dy = y(2)-y(1);
c = 0.2;
%dy = log(30)/(N-1);
[call, put] = blsprice(S,X,r,T,sigma);
for j = 1:N
    y(j) = (j-1)*dy;
end
phi = zeros(N,N);
for i = 1:N
    for j = 1:N
        phi(i,j) = sqrt((y(i) - y(j))^2 + c^2);
    end
end
L =  zeros(N,N);
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = (y(i)-y(j))/(sqrt((y(i) - y(j))^2 + c^2));
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 1/sqrt((y(i) - y(j))^2 + c^2)-(y(i)-y(j))^2/(sqrt((y(i) - y(j))^2 + c^2))^3;
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = max(X-S,0);
%alphao = zeros(N,1);
alphao = L\Uo';
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
    for j=2:N
        Unn(j,n)=Un(j);
    end
    alpha = L\Un;
    alphao = alpha;
end
figure(1)
plot(S,Unn(:,1),'-g',S,Unn(:,round(M/2)),'-b')
axis([0 18 0 9])
figure(2)
mesh(t,S,Unn)


% plot(S,Un,'k',S,put,'o')
% grid on
% t = dt:dt:T;
% for j = 1:M
%     F1 = -dt*P*alphao;
%     F2 = -dt*P*(alphao + 0.5*F1);
%     alpha = alphao + 0.5*(F1 + F2);
% %     Un(:,j) = L*alpha;
% %     Un(1) = X*exp(-r*j*dt);
% %     alpha = L\Un;
%     alphao = alpha;
% end
% plot(S,Un,'k',S,put,'o')
% grid on











