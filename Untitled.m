clear;
clc;
% ‰»Îx=[];

Y=x(:,1)
X=x(:,2)
g=[0.1];
n=10;
d=1;
t=[1:n]';
X1=cumsum(X);
Y1=cumsum(Y);
Z1=0.5*Y1(1:end-1)+0.5*Y1(2:end);
b=zeros(n-1,d);
for i=1:(n-1)
    for j=1:d
        for k=1:(i+1)
            b(i,j)=b(i,j)+g(j)^(i+1-k)*X1(k,j);
        end
    end
end
y=Y(2:end);
c=[2:n]'-1/2;
o=ones(n-1,1);
B=[-Z1,b,c,o];
B=B(1:(end),:);
pa=(B'*B)\B'*y;
pa2=pa(2:end-2);
mu1=1/(1+0.5*pa(1));
mu2=(1-0.5*pa(1))/(1+0.5*pa(1));
mu3=pa(end-1)/(1+0.5*pa(1));
mu4=(pa(end)-0.5*pa(end-1))/(1+0.5*pa(1));
Y1f(1)=Y(1);
for k=2:n
    Y1f(k,1)=mu1*b(k-1,:)*pa2+mu2*Y1f(k-1)+mu3*k+mu4;
end
Yf(1)=Y(1);
Yf(2:n,1)=diff(Y1f);
plot(t,Y,t,Yf)