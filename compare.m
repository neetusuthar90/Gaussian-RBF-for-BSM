clc;close all; clear all;
X = 100; r = 0.1; sigma = 0.30;
T = 1;
N =101;
M = 100;
y = linspace(0,6,N);
S = linspace(1,exp(6),N);
dt = T/M;
dy = y(2)-y(1);
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
lo = max(X-S,0);
bo = pinv(L)*lo';
b = zeros(N,1);
lnn = zeros(N,M);
for l = 1:M
    F1 = -dt*P*bo;
    F2 = -dt*P*(bo + 0.5*F1);
    F3 = -dt*P*(bo + 0.5*F2);
    F4 = -dt*P*(bo + F3);
    b = bo + (F1 + 2*F2 + 2*F3 + F4)/6;
    ln = L*b;
    lnn(1,l) = X*exp(-r*l*dt);
    for j=2:N
        lnn(j,l)=ln(j);
    end
    b = pinv(L)*ln;
    bo = b;
end
% A = [80;85;90;95;100;105;110;115;120;140;160;180];
% bino = [20.2689; 16.3467; 13.1228; 10.4847; 8.3348; 6.6071; 5.2091; 4.0976; 3.2059; 1.1789; 0.4231; 0.1501; 0.0529];
% plot(A,bino,'-b','')

% % plot(S,put,'-b','linewidth',1)
% % hold on
% % legend('European option')
% % xlabel('Stock price')
% % ylabel('Option value')
% plot (S,Unn(:,M),'-r','linewidth',1)
% hold on
% % legend('European option','American option')
% % plot(S,Uo,'--k','linewidth',1)
% % axis([0 180 0 100])
% % hold off
% % legend('European option','American option','Option''s payoff')
% plot(S,lnn(:,M),':g','linewidth',1)
% grid on
% % hold on
% % legend('Option''s payoff','put option', 'European option')

figure(2)
surf(lnn)

% xlabel('Time')
% ylabel('Stock price')
% zlabel('Option value')

