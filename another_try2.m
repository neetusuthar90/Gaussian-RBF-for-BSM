clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
theta = 0.5;
T = 0.5;
N = 141;
M = 30;
S = linspace(1,30,N);
dt = T/M;
y = log(S);
dy = y(2)-y(1);
%dy = 0.024;
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
%P=r*(eye(N))-(0.5*(S.*S)*sigma^2).*(pinv(L)*Lyy)-(r*S).*(pinv(L)*Ly);
%Uo = zeros(N,1);
Uo = max(X-S,0);
%plot(S,Uo)
%To = zeros(N,1);
alphao = GAUSS_ELIM(phi,Uo);
%To = L\Uo';
alpha = zeros(N,1);
P1 = (eye(N)+(1-theta)*dt*P);
P2 = (eye(N)-theta*dt*P);
Un = zeros(N,1);
Unn = zeros(N,M);
t = (1:M)*dt;
for j = 1:M
    t = T-t;
    b = P1\P2;
    alpha = b*alphao;
    Un = phi*alpha;
    for i = 1:N
        Unn(i,j) = Un(i);
    end
    
    alphao = alpha;
end
figure(1)
plot(S,Unn(:,1),'-g',S,Unn(:,round(M/2)),'-b')
    
figure(2)
mesh(t,S,Unn)