clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5;
N = 121;
M = 30;
smin = 1;
smax = 30;
S = linspace(smin,smax,N);
dt = T/M;
y = log(S);
dy = y(2) - y(1);
c = 1.5;
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
% A = [7.7531;
%     5.7531;
%     3.7532;
%     1.7987 ;
%       0.9880;
%     0.4420 ;
%     0.1606;
%     0.0483 ;
%       0.0124;
%     0.0028];
% B = [7.7751;
%     5.7751;
%     3.7842;
%     1.8937;
%     0.9549;
%     0.2545;
%     0.0812;
%     0.0210;
%     0.0056;
%     0.0020];
% tmp = 0;
% for i = 1:10
%     sum = tmp + (abs(A(i)-B(i)))^2;
%     tmp = sum;
% end
% RMSE = sqrt(sum/10);
tmp = 0;
for i = 1:N
    sum = tmp + (abs(put(i)-Unn(i,M)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/N);

figure(1)
plot(S,put,'-b','linewidth',1)
hold on
grid on
legend('Exact Solution')
axis([1 16 0 9])
xlabel('Asset price')
ylabel('Option value')
plot(S,Uo,'--k','linewidth',1)
hold on
legend('Exact Solution','Option value at t = 0.0')
plot(S,Unn(:,round(M/2)),'-g','linewidth',1)
hold on
legend('Exact Solution','Option value at t = 0.0','Option value at t = 0.25')
plot(S,Unn(:,M),'-m','linewidth',1)
hold on
hold off
legend('Exact Solution','Option value at t = 0.0','Option value at t = 0.25','Option value at t = 0.5')

% figure(2)
% surf(Unn)
% colormap cool;
% xlim([0 0.5])



