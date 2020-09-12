clc;close all; clear all;
X = 100; r = 0.1; sigma = 0.30;
T = 1;
N =101;
M = 100;
y = linspace(0,6,N);
S = linspace(1,exp(6),N);
dt = T/M;
dy = y(2)-y(1);
c = 4*dy;
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
Un =zeros(N,1);
t = (1:M)*dt;
for n = 1:M
    F1 = -dt*P*alphao;
    F2 = -dt*P*(alphao + 0.5*F1);
    alpha = alphao + 0.5*(F1 + F2);
    Un = L*alpha;
    for i = 1:N
        Un(i) = max(X-S(i), Un(i));
    end
    for j = 1:N
        Unn(j,n) = Un(j);
    end
    alpha = L\Un;
    alphao = alpha;
end
figure(1)
plot(S,Unn(:,1),'-g',S,Unn(:,round(M/2)),'-k',S,Unn(:,M),'-r')

figure(2)
surf(Unn)















