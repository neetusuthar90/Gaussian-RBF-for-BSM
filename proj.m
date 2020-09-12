clc
clear all
N=7; %spatial
M=30; %temporel
T=0.5; % Expiry
sigma=0.2;
r=0.05;
X=10;
S=linspace(1,30,N);
yy=log(S);
dyy=yy(2)-yy(1);
tt=linspace(0,T,M);
dtt=tt(2)-tt(1);
[t,y]=meshgrid(tt,yy);
phi=zeros(N,N);
U=zeros(N,M);
U(:,1)=max(X-exp(yy),0);
 for i=1:N
     phi(i,:)=sqrt((yy-yy(i)).^2+(4*dyy)^2);
 end
L=phi;
Ly=zeros(N,N);
for i=1:N
    for j=1:N
       Ly(i,j)=(yy(i)-yy(j))/L(i,j);
    end
end
Lyy=zeros(N,N);
for i=1:N
    for j=1:N
       Lyy(i,j)=(1/L(i,j))-((yy(i)-yy(j))^2/(L(i,j))^3);
    end
end

P=r*(eye(N))-(0.5*sigma^2)*(L\Lyy)-(r-(0.5*sigma^2))*(L\Ly);
[E,DV]=eig(P);
K=E\(L\U(:,1));
alpha=zeros(N,M);
Sum=zeros(N,1);
    for j=1:M
        for i=1:N
            Sum=Sum+K(i)*exp(DV(i,i)*t(1,j))*E(:,i);
           
        end
         alpha(:,j)=Sum;
    end
alpha=real(alpha);
sum=0;
for j=2:M
    for i=1:N
        for k=1:N
        sum=sum+alpha(k,j)*L(i,k);
        end
        U(i,j)=sum;
    end
end
