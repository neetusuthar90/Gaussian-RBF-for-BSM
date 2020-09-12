clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20; theta = 0.5;
T = 0.5;
N = 121;
M = 30;
S = linspace(1,30,N);
y = log(S);
dt = T/M;
dy = y(2)-y(1);
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
L =  zeros(N,N);
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = -2*c^2*(y(i)-y(j))*exp(-(c*(y(i)-y(j)))^2);
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 2*c^2*exp(-(c*(y(i)-y(j)))^2)*(2*c^2*(y(i)-y(j))^2 - 1);
    end
end
P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = max(X-S,0);
%alphao = zeros(N,1);
alphao = pinv(L)*Uo';
alpha = zeros(N,1);
P1 = (eye(N)+(1-theta)*dt*P);
P2 = (eye(N)-theta*dt*P);
Un = zeros(N,1);
Unn = zeros(N,M);
t = (1:M)*dt;
b = pinv(P1)*P2;
for j = 1:M
    alpha = b*alphao;
    Un = L*alpha;
    %Unn(1,j) = X*exp(-r*j*dt);
    for i = 1:N
        Unn(i,j) = Un(i);
    end
    alphao = alpha;
end


tmp = 0;
for i = 1:N
    sum = tmp + (abs(put(i)-Un(i)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/N);

% A = [8.7531;
%     7.7531;
%     5.7531;
%     3.7532;
%     1.7987 ;
%     0.4420 ;
%     0.0483 ;
%     0.0028 ;
%     0.0001];
% B =  [8.7531;
%     7.7531 ;
%     5.7531 ;
%     3.7533 ;
%     1.7996;
%     0.4411;
%     0.0487;
%     0.0037;
%     0.0008];
% for i = 1:9
%     sum = 0;
%     sum = sum + (abs(A(i)-B(i)))^2;
% end
% RMSE = sqrt(sum/9);
figure(1)
plot(S,put,'-b','linewidth',1)
hold on
grid on
legend('Exact solution')
axis([1 16 0 9])
xlabel('Stock price S')
ylabel('Option value')
plot(S,Unn(:,1),'--g','linewidth',1)
hold on
legend('Exact solution','Option value at t = 0')
plot(S,Unn(:,round(M/2)),'--m','linewidth',1)
hold on
legend('Exact solution','Option value at t = 0','Option value at t = 0.25')
plot(S,Unn(:,M),'-r','linewidth',1)
hold off
legend('Exact solution','Option value at t = 0','Option value at t = 0.25','Option value at t = 0.5')

% % 
% 
% figure(2)
% mesh(t,S,Unn)