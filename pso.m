function [ g , gb ] = pso(X,Y,N,d,T,c1,c2,w,xmax,xmin,vmax,vmin )
% function [ g , gb ] = pso(X,Y,s,d,T,c1,c2,w,xmax,xmin,vmax,vmin )
% TDAGM优化问题的粒子群算法
% X为相关因素数据序列
% Y为系统特征序列
% 输入：X,Y,s,d,T,c1,c2,w,xmax,xmin,vmax,vmin
% N,d,T分别为粒子群个数 自变量维数 循环次数
% c1,c2,w分别为学习参数1 学习参数2 惯性系数
% xmax,xmin为解空间边界
% vmax,vmin为速度限制
% 输出：
% g为参数估计值
% gb为目标函数进化值
%%
%%%%%%%%%%%%%%%%%初始化种群个体%%%%%%%%%%%%%%%%
xdiag=diag(xmax-xmin);
xminn=repmat(xmin,1,N);
x=xdiag*rand(d,N)+xminn;
vdiag=diag(vmax-vmin);
vminn=repmat(vmin,1,N);
v=vdiag*rand(d,N)+vminn;
%%
%%%%%%%%%%%初始化个体最优位置和最优值%%%%%%%%%%%
p=x;
for k=1:N
    pb(k)=fitc(X,Y,x(:,k));
end
%%%%%%%%%%初始化全局最优位置和最优值%%%%%%%%%%%%
[gb(1),gi]=min(pb);
g=x(:,gi);
%%
%%%%%%%迭代直至满足精度或达到最大迭代次数%%%%%%%%
for i = 2:T
    %%
        %%%%%%%%%%%%更新位置和速度%%%%%%%%%%%%%%
        gn=repmat(g,1,N);
        r=rand(N,1);
        r1=diag(r);
        r=rand(N,1);
        r2=diag(r);
        v=w*v+c1*(p-x)*r1+c2*(gn-x)*r2;
        x=x+v;
    %%
        %%%%%%%%5%%%%边界条件处理%%%%%%%%%%%%%%%
        for j=1:d
        flag=x(j,:)>xmax(j)|x(j,:)<xmin(j);
        l=sum(flag(:));
        x(j,flag)=rand(l,1)*(xmax(j)-xmin(j))+xmin(j);
        flag=v(j,:)>vmax(j)|v(j,:)<vmin(j);
        l=sum(flag(:));
        v(j,flag)=rand(l,1)*(vmax(j)-vmin(j))+vmin(j);
        end
    %%
        %%%%%%%更新个体最优位置和最优值%%%%%%%%%%
        for k=1:N
        pm(k)=fitc(X,Y,x(:,k));
        end
        flag=pm<pb;
        pb(flag)=pm(flag);
        p(:,flag)=x(:,flag);
        %%%%%%%更新全局最优位置和最优值%%%%%%%%%%
        [gb(i),gi]=min(pb);
        g=p(:,gi);
end
end