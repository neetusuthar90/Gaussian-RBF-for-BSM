clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20;
T = 0.5;
N = 160;
M = 1600;
Smin = 0;
Smax = 30;
dS = (Smax-Smin)/N;
dt = T/M;
V(1:N+1,1:M+1)= 0;
S = Smin+(0:N)*dS;
tau = (0:M)*dt;
V(1,:) = Smax - X*exp(-r*tau);
V(N+1,:) = 0;
V(:,M+1)= max(X-S,0);
for j = 1:M
    for n = 2:N
        V(n,j+1)=0.5*dt*(sigma*sigma*n*n-r*n)*V(n-1,j) + (1 - dt*(sigma*sigma*n*n+r))*V(n,j)+0.5*dt*(sigma*sigma*n*n+r*n)*V(n+1,j);
    end
end
figure(1)
plot(S,V(:,1),'r-',S,V(:,round(M/2)),'g-',S,V(:,M+1),'b-')
xlabel('S')
ylabel('V')

figure(2)
mesh(tau,S,V)
