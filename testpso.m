%% ��������
% clc;
% load foodstuff.mat;
s=7;
XC=X(1:s,:);
YC=Y(1:s);
%% ����Ⱥ�Ż�����
N=100;
d=3;
T=200;
c1=1.5;
c2=1.5;
w=0.8;
xmax=[1;1;1];
xmin=[0;0;0];
vmax=[0.1;0.1;0.1];
vmin=[-0.1;-0.1;-0.1];
%%
%%%%%%%%%%%%%%%��ʼ����Ⱥ����%%%%%%%%%%%%%%%%
xdiag=diag(xmax-xmin);
xminn=repmat(xmin,1,N);
x=xdiag*rand(d,N)+xminn;
vdiag=diag(vmax-vmin);
vminn=repmat(vmin,1,N);
v=vdiag*rand(d,N)+vminn;
%%
%%%%%%%%%��ʼ����������λ�ú�����ֵ%%%%%%%%%%%
p=x;
for k=1:N
    pb(k)=fitc(XC,YC,x(:,k));
end
%%%%%%%%%%��ʼ��ȫ������λ�ú�����ֵ%%%%%%%%%%%%
[gb(1),gi]=min(pb);
g=x(:,gi);
%%
%%%%%%%����ֱ�����㾫�Ȼ�ﵽ����������%%%%%%%%
for i = 2:T
    %%
        %%%%%%%%%%%%����λ�ú��ٶ�%%%%%%%%%%%%%%
        gn=repmat(g,1,N);
        r=rand(N,1);
        r1=diag(r);
        r=rand(N,1);
        r2=diag(r);
        v=w*v+c1*(p-x)*r1+c2*(gn-x)*r2;
        x=x+v;
    %%
        %%%%%%%%5%%%%�߽���������%%%%%%%%%%%%%%%
        for j=1:d
        flag=x(j,:)>xmax(j)|x(j,:)<xmin(j);
        l=sum(flag(:));
        x(j,flag)=rand(l,1)*(xmax(j)-xmin(j))+xmin(j);
        flag=v(j,:)>vmax(j)|v(j,:)<vmin(j);
        l=sum(flag(:));
        v(j,flag)=rand(l,1)*(vmax(j)-vmin(j))+vmin(j);
        end
    %%
        %%%%%%%���¸�������λ�ú�����ֵ%%%%%%%%%%
        for k=1:N
        pm(k)=fitc(XC,YC,x(:,k));
        end
        flag=pm<pb;
        pb(flag)=pm(flag);
        p(:,flag)=x(:,flag);
        %%%%%%%����ȫ������λ�ú�����ֵ%%%%%%%%%%
        [gb(i),gi]=min(pb);
        g=p(:,gi);
end
%%���
[n,d]=size(X);
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
y=Y(2:s);
c=[2:n]'-1/2;
o=ones(n-1,1);
B=[-Z1,b,c,o];
B=B(1:(s-1),:);
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