clc;
clear all;
close all;
X = 10;
r = 0.05;
N = 21;
M = 30;
sigma = 0.20;
t = 0;
T = 0.5; %year
Smin = 1;
Smax = 30;
y = zeros(1,N);
S = linspace(Smin,Smax,N);
t = linspace(0,T,M);
[t,y]=meshgrid(t,y);
Y = log(S);
dy = log(Smax)/(N-1);
dt = T/M;
[call,Put] = blsprice(S,X,r,T,sigma);
for j = 1:N
    y(j) = (j-1)*dy;
end
for i = 2:N
    d(i) = y(i) - y(i-1);
end
c = 0.5;
for k = 2:N
    if d(k) < c
      c = d(k);
    end
end
phi = zeros(N,N);
for i=1:N
     phi(i,:)=sqrt((Y-y(i)).^2+c^2);
 end
L = zeros(N,N);
for i = 1:N
    for j = 1:N
        L(i,j) = sqrt((y(i)-y(j))^2 + c^2);
    end
end
Ly = zeros(N,N);
for i = 1:N
    for j = 1:N
        Ly(i,j) = (y(i)-y(j))/sqrt((y(i)-y(j))^2 + c^2);
    end
end
Lyy = zeros(N,N);
for i = 1:N
    for j = 1:N
        Lyy(i,j) = 1/sqrt((y(i)-y(j))^2 + c^2) - (y(i)-y(j))^2/((y(i)-y(j))^2 + c^2)^(1.5);
    end
end
P = r*(eye(N)) - (0.5*sigma^2)*(L\Lyy) - (r-0.5*sigma^2)*(L\Ly);
% %Represent eigenvalue of P by l
U=zeros(N,M);
U(:,1)=max(X-exp(Y),0);
[W,D] = eig(P);
K=W\(L\U(:,1));
% Analytical time integration formula
alpha=zeros(N,M);
Sum=zeros(N,1);
    for j=1:M
        for i=1:N
            Sum=Sum+K(i)*exp(D(i,i)*t(1,j))*W(:,i); 
        end
         alpha(:,j)=Sum;
    end
alpha=real(alpha);
sum=0;
for j=2:M
    for i=1:N
        for k=1:N
        sum=sum+alpha(k,j)*phi(i,k);
        end
        U(i,j)=sum;
    end
end




% Explicit second order backward time integration scheme
% here A is initial vector alpha(o)
% alpha0 = zeros(M,M);
% alpha = zeros(M,M);
% alpha0 = L\U(:,1);
% Un = zeros(M,M);
% for n = 2:M
%     F1 = -dt*(P*alpha0);
%     F2 = -dt*P*(alpha0+ 0.5*F1);  
%     alpha = alpha0 + 0.5*(F1 + F2);
%     Un = L*alpha;
%     Un(1) = X*exp(-r*n*dt);
%     alpha = L\Un;
%     alpha0 = alpha;
% end 
% display(Un)

% theta = 0.5;
% P1 = (eye(N) + (1-theta)*dt*P);
% P2 = (eye(N) - theta*dt*P);
% alpha0 = zeros(N,1);
% an = zeros(N,1);
% alpha0 = L\U0;
% for j = 1:M
%     an = P1\P2*alpha0;
%     alpha0 = an;
% end
% an
