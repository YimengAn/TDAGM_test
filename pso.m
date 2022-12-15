function [ g , gb ] = pso(X,Y,N,d,T,c1,c2,w,xmax,xmin,vmax,vmin )
% function [ g , gb ] = pso(X,Y,s,d,T,c1,c2,w,xmax,xmin,vmax,vmin )
% TDAGM�Ż����������Ⱥ�㷨
% XΪ���������������
% YΪϵͳ��������
% ���룺X,Y,s,d,T,c1,c2,w,xmax,xmin,vmax,vmin
% N,d,T�ֱ�Ϊ����Ⱥ���� �Ա���ά�� ѭ������
% c1,c2,w�ֱ�Ϊѧϰ����1 ѧϰ����2 ����ϵ��
% xmax,xminΪ��ռ�߽�
% vmax,vminΪ�ٶ�����
% �����
% gΪ��������ֵ
% gbΪĿ�꺯������ֵ
%%
%%%%%%%%%%%%%%%%%��ʼ����Ⱥ����%%%%%%%%%%%%%%%%
xdiag=diag(xmax-xmin);
xminn=repmat(xmin,1,N);
x=xdiag*rand(d,N)+xminn;
vdiag=diag(vmax-vmin);
vminn=repmat(vmin,1,N);
v=vdiag*rand(d,N)+vminn;
%%
%%%%%%%%%%%��ʼ����������λ�ú�����ֵ%%%%%%%%%%%
p=x;
for k=1:N
    pb(k)=fitc(X,Y,x(:,k));
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
        pm(k)=fitc(X,Y,x(:,k));
        end
        flag=pm<pb;
        pb(flag)=pm(flag);
        p(:,flag)=x(:,flag);
        %%%%%%%����ȫ������λ�ú�����ֵ%%%%%%%%%%
        [gb(i),gi]=min(pb);
        g=p(:,gi);
end
end