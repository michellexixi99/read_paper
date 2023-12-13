%Generalizing RNA velocity to transient cell states through dynamical modeling

clc;clear all;close all;

%normalized data
% s_obs=load('spliced_data.txt');
% u_obs=load('unspliced_data.txt');


%initialize at ti
for ch=1:10%5
s1=load('s1.mat'); 
s2=load('s2.mat'); 
s1= s1.s1;
s2= s2.s2;
s_obs=[s1,s2];
u_obs=load('u_obs.mat');
% s_obs=s_obs.s_obs;
u_obs=u_obs.u_obs;

lb=ch;ub=100-ch;
beta=1;
% lb=5;ub=95;
s_obs=s_obs(:,2);
u_obs=u_obs(:,2);
indi1=find(s_obs>0);
s_obs=s_obs(indi1);
u_obs=u_obs(indi1);
indi2=find(u_obs>0);
s_obs=s_obs(indi2);
u_obs=u_obs(indi2);

gamap=lf(s_obs,u_obs,lb,ub);  %
gamma=gamap*beta;
utk1=[];
stk1=[];
std_u=1;
std_s=1;

u_pre=ones(length(u_obs),1);
s_pre=ones(length(u_obs),1);
u00=ones(2,1)*0.9;
s00=ones(2,1)*0.9;

tao=ones(length(u_obs),1)*0.1;%==?===========================================

t000=[1,1];
% s_obs=s_obs%/max(s_obs);
% u_obs=u_obs%/max(u_obs);
%innitialize aplha
alpha2k(1)=min(max(s_obs),1);
alpha2k(2)=0;


for i=1:length(s_obs)
    if u_obs(i)-gamap*s_obs(i)>0
        k(i)=1;
        alpha(i)=real(alpha2k(1));
        u0(i)=u00(1);
        s0(i)=s00(1);
    else
        k(i)=0;
        alpha(i)=real(alpha2k(2));
        u0(i)=u00(2);
        s0(i)=s00(2);
    end
        tao(i)=real(solve_tao(alpha(i),gamma,beta,u_obs(i),s_obs(i),s0(i),u0(i)));
        [u_pre(i),s_pre(i)]=eq4(alpha(i), beta, gamma, tao(i), u0(i), s0(i));
end
diff=1;
inte=1;
max_iter=100;    
t=zeros(length(u_obs),1);
while inte<=max_iter
    t00=t000;
% while    diff>1e-2
%     alpha00=alpha2k;
    k1=find(k==1);
    k0=find(k==0);
    %update time

    t(k1)=t00(1)+tao(k1);
    t(k0)=t00(2)+tao(k0);

    beta=solveforbata(u_obs,u_pre,s_obs,s_pre,u0,alpha',gamma,tao,s0);
    
    alpha2k(1)=solveforalpha(u_obs(k1),u_pre(k1),s_obs(k1),s_pre(k1),u0(k1),gamma,t(k1), beta,tao(k1), s0(k1)); 
    t00(1)=solvefort0(u_obs(k1),u_pre(k1),s_obs(k1),s_pre(k1),u0(k1),alpha2k(1),beta,s0(k1),tao(k1),gamma,t(k1));
    alpha2k(2)=solveforalpha(u_obs(k0),u_pre(k0),s_obs(k0),s_pre(k0),u0(k0),gamma,t(k0), beta,tao(k0), s0(k0));
    t00(2)=solvefort0(u_obs(k0),u_pre(k0),s_obs(k0),s_pre(k0),u0(k0),alpha2k(2),beta,s0(k0),tao(k0),gamma,t(k0));
    gamma=solveforgamma(u0,alpha',beta,s0,tao,s_obs);
    gamap=gamma/beta;
%     gama=gamap;
    dft=t00-t000;
    [u00(1),s00(1)]=eq4(alpha2k(1), beta, gamma, abs(dft(1)), u00(1), s00(1));
    [u00(2),s00(2)]=eq4(alpha2k(2), beta, gamma, abs(dft(2)), u00(2), s00(2));

    
    alpha(k1)=alpha2k(1);
    alpha(k0)=alpha2k(2);
    u0(k1)=u00(1);
    u0(k0)=u00(2);
    s0(k1)=s00(1);
    s0(k0)=s00(2);

    %update k t
    for i=1:length(s_obs)
%        tao(i)=real(solve_tao(alpha(i),gamma,beta,u_obs(i),s_obs(i),s0(i),u0(i)));
       tao(i)=real(solve_tao(alpha(i),gamma,beta,u_obs(i),s_obs(i),s0(i),u0(i)));
       [u_pre1(i),s_pre1(i)]=eq4(alpha2k(1), beta, gamma, tao(i), u00(1), s00(1));
       [u_pre0(i),s_pre0(i)]=eq4(alpha2k(2), beta, gamma, tao(i), u00(2), s00(2));
       dis1(i)= (cal_dis(u_obs(i),u_pre1(i),std_u))^2+(cal_dis(u_obs(i),u_pre1(i),std_s))^2;
       dis2(i)= (cal_dis(u_obs(i),u_pre0(i),std_u))^2+(cal_dis(u_obs(i),u_pre0(i),std_s))^2;
       if dis1(i)>dis2(i)    
           k(i)=0;
       else
           k(i)=1;
       end
       [u_pre(i),s_pre(i)]=eq4(alpha(i), beta, gamma, tao(i), u0(i), s0(i));

    end
    

% diff=norm(alpha2k-alpha00,2);

% norm([u_pre-u_obs,s_pre-s_obs],2)
    xx_pre=[];
    xx_obs=[];
    xx_pre=[u_pre,s_pre];
    xx_obs=[u_obs,s_obs];
    corr_matrix = corrcoef(xx_pre, xx_obs);
    r = corr_matrix(1, 2);
    r2(inte) = r^2;
    inte=1+inte; 
end
r22(ch)=r2(inte-1);
alpha2k_all(ch,1:2)=alpha2k;
alpha2k=[];
beta_all(ch,1)=beta;
beta=[];
gamma_all(ch,1)=gamma;
gamma=[];
end

%alpha
plot(1:ch,alpha2k_all(1:ch,1));hold on;
plot(1:ch,alpha2k_all(1:ch,2));
plot([1 ch],[0.38913 0.38913])
yticks(gca,[0:0.2:1]);
ylim([0, 0.6]);
xlabel('Extreme quantiles','fontname','Times New Roman');
ylabel('Alpha','fontname','Times New Roman');
grid on;
box on;
legend('alpha{k=1}','alpha{k=0}','scVelo','fontname','Times New Roman');

%beta
plot(1:ch,beta_all(1:ch,1));hold on;
plot([1 ch],[2.981981 2.981981],'color',[0.93,0.69,0.13]);
yticks(gca,[0:1:4]);
xlabel('Extreme quantiles','fontname','Times New Roman');
ylabel('Beta','fontname','Times New Roman');
grid on;
box on;
legend('Beta','scVelo','fontname','Times New Roman');


%gamma
plot(1:ch,gamma_all(1:ch,1));hold on;
plot([1 ch],[0.26032 0.26032],'color',[0.93,0.69,0.13]);
yticks(gca,[0:0.2:1]);
ylim([0, 1]);
xlabel('Extreme quantiles','fontname','Times New Roman');
ylabel('Gamma','fontname','Times New Roman');
grid on;
box on;
legend('Gamma','scVelo','fontname','Times New Roman');


%r2
plot(1:ch,r22(1:ch));hold on;
yticks(gca,[0:0.2:1]);
ylim([0, 1]);
xlabel('Extreme quantiles','fontname','Times New Roman');
ylabel('R^2','fontname','Times New Roman');
grid on;
box on;
