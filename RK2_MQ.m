clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5;
N =81;
M = 30;
smin = 1;
smax = 30;
S = linspace(smin,smax,N);
ymin = log(smin);
ymax = log(smax);
dt = T/M;
y = linspace(ymin,ymax,N);
%dy = y(2)-y(1);
dy = log(30)/(N-1);
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
        Lyy(i,j) = 1/sqrt((y(i) - y(j))^2 + (4*dy)^2)-(y(i)-y(j))^2/(sqrt((y(i) - y(j))^2 + (4*dy)^2))^3;
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = max(X-S,0);
alphao = pinv(L)*(Uo');
alpha = zeros(N,1);
Unn = zeros(N,M);
for n = 1:M
    t = n*dt;
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    alpha = alphao + 0.5*(F1 + F2);
    Un = L*alpha;
    Unn(1,n) = X*exp(-r*t);
    for j=2:N
        Unn(j,n)=Un(j);
    end
    alpha = L\Un;
    alphao = alpha;
end
tmp = 0;
for i = 1:N
    sum = tmp + (abs(put(i)- Un(i)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/N);

figure(1)
plot(S,put,'-g','linewidth',1)
hold on
grid on
legend('Analytical option value')
plot(S,Unn(:,1),':r','linewidth',1)
hold on
legend('Analytical option value','Option value at t = 0')
plot(S,Unn(:,round(M/2)),'--r','linewidth',1)
hold on
legend('Analytical option value','Option value at t = 0','Option value at t = 0.25')
plot(S,Unn(:,M),'-r','linewidth',1)
hold on
legend('Analytical option value','Option value at t = 0','Option value at t = 0.25','Option value at t = 0.5')
axis([2 18 0 9])
hold off
xlabel('Stock price S')
ylabel('Option value')







