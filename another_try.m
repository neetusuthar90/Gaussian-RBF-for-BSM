clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5; theta = 0.5;
N = 61;
M = 5;
S = linspace(1,30,N);
dt = T/M;
y = linspace(log(1),log(30),N);
%dy = log(30)/(N-1);
dy = y(2)-y(1);
c = 1;
for j = 1:N
    y(j) = (j-1)*dy;
end
phi = zeros(N,N);
for i = 1:N
    for j = 1:N
        phi(i,j) = sqrt(abs(y(i) - y(j))^2 + 4*dy^2);
    end
end
L =  zeros(N,N);
Ly =  zeros(N,N);
Lyy =  zeros(N,N);
L = phi;
for i = 1:N
    for j = 1:N
        Ly(i,j) = (y(i)-y(j))/(sqrt(abs(y(i) - y(j))^2 + 4*dy^2));
    end
end
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 1/sqrt(abs(y(i) - y(j))^2 + 4*dy^2)-(y(i)-y(j))^2/(sqrt(abs(y(i) - y(j))^2 + 4*dy^2))^3;
    end
end
 P=r*(eye(N))-(0.5*sigma^2)*(pinv(L)*Lyy)-(r-(0.5*sigma^2))*(pinv(L)*Ly);
Uo = zeros(N,1);
Uo = max(X-exp(y),0);
%plot(S,Uo)
To = zeros(N,1);
To = GAUSS_ELIM(phi,Uo);
%To = L\Uo';
Tn = zeros(N,1);
P1 = zeros(N,N);
P2 = zeros(N,N);
P1 = (eye(N)+(1-theta)*dt*P);
P2 = (eye(N)-theta*dt*P);
Un = zeros(N,1);
for i = dt:dt:0.5
    b = inv(P1)*P2;
    Tn = b*To;
    Un = phi*To;
    To = Tn;
%     plot(S,Un)
%     grid on
end

    
    
    