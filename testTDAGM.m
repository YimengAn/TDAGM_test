[n,N]=size(X);
X1=cumsum(X);
Y1=cumsum(Y);
Z1=0.5*Y1(1:end-1)+0.5*Y1(2:end);
b=zeros(n-1,N);
for i=1:(n-1)
    for j=1:N
        for k=1:(i+1)
            b(i,j)=b(i,j)+l(j)^(i+1-k)*X1(k,j);
        end
    end
end
y=Y(2:end);
c=[2:n]'-1/2;
o=ones(n-1,1);
B=[-Z1,b,c,o];
p=(B'*B)\B'*y;
pb=p(2:end-2);
mu1=1/(1+0.5*p(1));
mu2=(1-0.5*p(1))/(1+0.5*p(1));
mu3=p(end-1)/(1+0.5*p(1));
mu4=0.5*p(end)/(1+0.5*p(1));
Yf(1)=Y(1);
for k=2:n
    Yf(k)=mu1*b(k-1,:)*pb+mu2*Yf(k-1)+mu3*(k-1/2)+mu4;
end
YY=Yf'