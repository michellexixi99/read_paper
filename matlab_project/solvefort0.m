function t00=solvefort0(u_obs,u_pre,s_obs,s_pre,u0,alpha,bata,s0,ttao,gamma,t)
syms t0

% ut=u0*exp(-bata*tao)+alpha/bata*(1-exp(-bata*tao));
% st=s0*exp(-gamma*tao)+alpha/gamma*(1-exp(-gamma*tao))+(alpha-bata*u0)/(gamma-bata)*(exp(-gamma*tao)-exp(-bata*tao));
% ua=diff(ut,alpha)
% sa=diff(st,alpha)
% ub=diff(ut,bata)
% sb=diff(st,bata)
% ug=diff(ut,gamma)
% sg=diff(st,gamma)
% ut0=diff(ut,t0)
% st0=diff(st,t0)


t0_xuanz=0.1:0.5:40;
for j=1:length(t0_xuanz)
    t0=t0_xuanz(j);
    tao=ttao+(t-ttao)-t0;
    for i=1:size(u_obs,1)
        [u_pre(i),s_pre(i)]=eq4(alpha, bata, gamma, tao(i), u0(i), s0(i));
    end
    
ut0 =bata.*u0'.*exp(-bata.*(tao)) - alpha.*exp(-bata.*(tao));
st0 =gamma.*s0'.*exp(-gamma.*(tao)) - alpha.*exp(-gamma.*(tao)) + ((bata.*exp(-bata.*(tao)) - gamma.*exp(-gamma.*(tao))).*(alpha - bata.*u0'))./(bata - gamma);
f(j)=sum(abs(ut0.*sum(u_obs-u_pre)+st0.*sum(s_obs-s_pre)));
end
[val,indi]=min(f);
t00=t0_xuanz(indi);

end