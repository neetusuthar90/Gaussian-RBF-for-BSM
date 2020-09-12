clc;
close all;
clear all;
X = 100; r = 0.1; sigma = 0.30;
T = 1;
N =101;
M = 100;
y = linspace(0,6,N);
S = exp(y);
%S = linspace(1,exp(6),N);
dt = T/M;
dy = y(2)-y(1);
%dy = log(30)/(N-1);
c = 6;
[call, put] = blsprice(S,X,r,T,sigma);
% bino = [20.2689; 16.3467; 13.1228; 10.4847; 8.3348; 6.6071; 5.2091; 4.0976; 3.2059; 1.1789; 0.4231; 0.1501; 0.0529];
phi = zeros(N,N);
for j = 1:N
    y(j) = (j-1)*dy;
end
for i = 1:N
    for j = 1:N
        phi(i,j) = exp(-c^2*(y(i)-y(j))^2);
    end
end
%L =  zeros(N,N);
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
t = (1:M)*dt;
Unn = zeros(N,M);
for n = 1:M
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    alpha = alphao + 0.5*(F1 + F2);
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
B = [20.2814;
  16.3552 ;
  13.1124 ;
 10.5051 ;
  8.3855;
  6.6652 ;
 5.2627;
  4.1108; 
  3.2419 ;
   1.2024;
   0.437 ;
  0.1641
  ];
tmp = 0;
for i = 1:12
    sum = tmp + (abs(A(i)-B(i)))^2;
    tmp = sum;
end
RMSE = sqrt(sum/12);
S1 = [80;85;90;95;100;105;110;115;120;140;160;180];
figure(1)
plot(S1,A,'--k','linewidth',1)
hold on
grid on
plot(S,Unn(:,M/4),':g','linewidth',1)
hold on
plot(S,Unn(:,round(M/2)),'-b','linewidth',1)
hold on
plot(S,Unn(:,M),'--r','linewidth',1)
hold on
axis([0 180 0 100])
hold off
xlabel('Stock price S')
ylabel('Option value')

legend('Binomial value','Option value at t = 0.25','Option value at t = 0.50','Option value at t = 1.0')

figure(2)
mesh(t,S,Unn)