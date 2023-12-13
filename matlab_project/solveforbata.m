function bataop=solveforbata(u_obs,u_pre,s_obs,s_pre,u0,alpha,gamma,tao,s0)
% syms bata

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
1

bata_xuanz=0.01:0.1:4;
for j=1:length(bata_xuanz)
    bata=bata_xuanz(j);
    for i=1:size(u_obs,1)
        [u_pre(i),s_pre(i)]=eq4(alpha(i), bata, gamma, tao(i), u0(i), s0(i));
    end

    ub=(alpha.*(exp(-bata.*(tao)) - 1))./bata.^2 - u0'.*exp(-bata.*(tao)).*(tao) + (alpha.*exp(-bata.*(tao)).*(tao))./bata;
    sb=-(u0'.*(exp(-bata.*(tao)) - exp(-gamma.*(tao))))./(bata - gamma) - ((exp(-bata.*(tao)) - exp(-gamma.*(tao))).*(alpha - bata*u0'))./(bata - gamma).^2 - (exp(-bata.*(tao)).*(tao).*(alpha - bata*u0'))./(bata - gamma);
    f(j)=abs(sum(ub*sum(u_obs-u_pre)+sb*sum(s_obs-s_pre)));
end
[val,indi]=min(f);
bataop=bata_xuanz(indi);

end