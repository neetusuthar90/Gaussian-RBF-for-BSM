clc;
close all;
clear all;
X = 100; r = 0.1; sigma = 0.3;
T = 1;
N = 121;
M = 100;
y = linspace(0,6,N);
S = exp(y);
%S = linspace(1,exp(6),N);
dt = T/M;
dy = y(2)-y(1);
%dy = log(30)/(N-1);
c = 7;
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
alphao = pinv(L)*Uo';
alpha = zeros(N,1);
Unn = zeros(N,M);
t = (1:M)*dt;
for n = 1:M
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    F3 = -dt*P*(alphao + 0.5*F2);
    F4 = -dt*P*(alphao + F3);
    alpha = alphao + (F1 + 2*F2 + 2*F3 + F4)/6;
    Un = L*alpha;
    for i = 1:N
        Un(i) = max(X- S(i), Un(i));
    end
    for j = 1:N
        Unn(j,n) = Un(j);
    end
    alpha = L\Un;
    alphao = alpha;
end
A = [20.2689; 
   16.3467 ;
  13.1228;
  10.4847;
  8.3348 ;
   6.6071 ;
   5.2091 ;
   4.0976 ;
   3.2059 ;
   1.1789 ;
   0.4231 ;
   0.1502];
B = [ 20.2631;
    16.3341;
    13.0963;
    10.4911;
    8.375;
    6.6560;
    5.2717;
    4.1012;
    3.2341;
    1.2003;
    0.4376;
    0.1651 ];
tmp = 0;
for i = 1:12
    sum = tmp + (abs(A(i)-B(i)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/12);

S1 = [80;85;90;95;100;105;110;115;120;140;160;180];
figure(1)
plot(S,Unn(:,M/4),':g','linewidth',1)
hold on
plot(S,Unn(:,round(M/2)),'--k','linewidth',1)
hold on
plot(S,Unn(:,M),'-r','linewidth',1)
hold on
% plot(S1,A,'--y','linewidth',1)
% hold on
axis([0 180 0 100])
hold on
xlabel('Stock price S')
ylabel('Option value')
legend('Option value at t = 0.25','Option value at t = 0.50','Option value at t = 1.0','Binomial value')
plot(S,Uo,'--b','linewidth',1)
hold off
legend('Option value at t = 0.25','Option value at t = 0.50','Option value at t = 1.0','Binomial method','Option Payoff')


% figure(2)
% mesh(t,S,Unn)



