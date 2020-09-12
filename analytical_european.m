clc;
clear all;
close all;
X = 10;
r = 0.05;
sigma = 0.20;
T = 0.5; %year
Smin = 1;
Smax = 30;
N = 21;
t = 0;
M = 30;
dy = log(Smax)/(N-1);
dt = T/M;
y = zeros(1,N);
for j = 1:N
    y(j) = (j-1)*dy;
   % fprintf('y(%d) = %f\n',j,y(j));
end
for i = 2:N
    d(i) = y(i) - y(i-1);
    %fprintf('d(%d) = %f\n',i,d(i));
end
c = 0.5;
for k = 2:N
    if d(k) < c
      c = d(k);
    end
end
display(c)
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
%let V be the inverse of L
V = gauss_elimination(L);
%Find the matrix P
I = eye(N);
P = r*I - (0.5*sigma^2)*(V*Lyy) - (r-0.5*sigma^2)*(V*Ly);
%display(P)
for i = 1:N
    U(i) = max(X-exp(y(i)),0);
end
%the initial vector U0 of U(yi,T)
U0 = U(1);

%Represent eigenvalue of P by l
l = eig(P);
[W,D] = eig(P);
%Inverse of matrix W is represented by Z
Z = gauss_elimination(W);
k = Z*V*U0;
% Analytical time integration formula
A = zeros(N,1);
for n = 1:M
    for i = 1:N
        x = W(:,i);
        y = k(:,i);
        A = A + exp(l(i)*t)*(y.*x);
        t = T-n*dt;
    end
end
Ar = real(A);
S = 2;
Y = log(S);
U = 0;
for j = 1:N
    f = sqrt((Y-y(j))^2 + c^2);
    U = U + Ar(j)*f;
end
display(U)