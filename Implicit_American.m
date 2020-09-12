clc;
close all;
clear all;
X = 100; r = 0.1; sigma = 0.30;
T = 1;
N =101;
M = 100;
theta = 0.5;
y = linspace(0,6,N);
S = linspace(1,exp(6),N);
%S = linspace(1,exp(6),N);
dt = T/M;
dy = y(2)-y(1);
%dy = log(30)/(N-1);
c = 6;
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
%plot(S,Uo)
%To = zeros(N,1);
%alphao = GAUSS_ELIM(phi,Uo);
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
    for i = 1:N
        Un(i) = max(X- S(i), Un(i));
    end
        for i = 1:N
            Unn(i,j) = Un(i);
        end
     alpha = L\Un;
    alphao = alpha;
end
B = [20.2890;
    16.3712;
    13.1291;
    10.4905;
    8.3471;
    6.6137;
    5.2210;
    4.1210;
    3.2187;
    1.1829;
    0.4651;
    0.1632];
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
tmp = 0;
for i = 1:12
    sum = tmp + (abs(A(i)-B(i)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/12);
S1 = [80;85;90;95;100;105;110;115;120;140;160;180];
figure(1)
plot(S,Uo,'--r','linewidth',1)
hold on
grid on
legend('Option''s payoff')
axis([0 180 0 100])
xlabel('Stock price S')
ylabel('Option value')
plot(S,Unn(:,M/4),'-b','linewidth',1)
hold on
legend('Option''s payoff','Option value at t = 0.25')
plot(S,Unn(:,round(M/2)),'-m','linewidth',1)
hold on
legend('Option''s payoff','Option value at t = 0.25','Option value at t = 0.50')
plot(S,Unn(:,M),'-k','linewidth',1)
hold on
legend('Option''s payoff','Option value at t = 0.25','Option value at t = 0.50','Option value at t = 1.0')
plot(S1,A,'--g','linewidth',1)
hold on
hold off
legend('Option''s payoff','Option value at t = 0.25','Option value at t = 0.50','Option value at t = 1.0','Binomial value');
