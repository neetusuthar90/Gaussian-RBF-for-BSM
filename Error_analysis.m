clc;close all; clear all;
X = 10; r = 0.05; sigma = 0.20; theta = 0.5;
T = 0.5;
N = 81;
M = 100;
S = linspace(1,30,N);
y = log(S);
dt = T/M;
dy = y(2)-y(1);
c = 1.53;
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
%uo = zeros(N,1);
uo = pinv(L)*Uo';
vo = uo;
wo = uo;
u = zeros(N,1);
v = u;
w = u;
Un = zeros(N,1);
Vn = zeros(N,1);
Wn = zeros(N,1);
Unn = zeros(N,M);
Vnn = zeros(N,M);
Wnn = zeros(N,M);
t = (1:M)*dt;
P1 = (eye(N)+(1-theta)*dt*P);
P2 = (eye(N)-theta*dt*P);
b = pinv(P1)*P2;
for n = 1:M
    
    %Implicit
    u = b*uo;
    Un = L*u;
    for i = 1:N
        Unn(i,n) = Un(i);
    end
    uo = u;
end
for n = 1:M
    %RK2
    f1 = -dt*P*vo;
    f2 = -dt*P*(vo + 0.5*f1);
    v = vo + 0.5*(f1 + f2);
    Vn = L*v;
    Vnn(1,n) = X*exp(-r*n*dt);
    for j=2:N
        Vnn(j,n)=Vn(j);
    end
    v = L\Vn;
    vo = v;
end   
for n = 1:M
    %RK4
    F1 = -dt*P*wo;
    F2 = -dt*P*(wo + 0.5*F1);
    F3 = -dt*P*(wo + 0.5*F2);
    F4 = -dt*P*(wo + F3);
    w = wo + (F1 + 2*F2 + 2*F3 + F4)/6;
    Wn = L*w;
    Wnn(1,n) = X*exp(-r*n*dt);
    for j=2:N
        Wnn(j,n)=Wn(j);
    end
    w = L\Wn;
    wo = w;
end

figure(1)
plot(S,Unn(:,1),'b',S,Unn(:,round(M/2)),'b',S,Unn(:,M),'b','linewidth',1)
hold on
plot(S,Vnn(:,1),'r',S,Vnn(:,round(M/2)),'r',S,Vnn(:,M),'r','linewidth',1)
hold on
plot(S,Wnn(:,1),'g',S,Wnn(:,round(M/2)),'g',S,Wnn(:,M),'g','linewidth',1)
hold on
axis([2 18 0 9])
hold off
xlabel('Stock price S')
ylabel('Option value')

