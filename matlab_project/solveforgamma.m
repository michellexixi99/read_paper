% function gamm=solveforgammaa(u_obs,u_pre,s_obs,s_pre,u0,tao,alpha,bata,s0,t,t0)
function gamm=solveforgamma(u0,alpha,bata,s0,tao,s_obs)

gammaa_xuanz=0.01:0.05:0.5;
f=zeros(length(gammaa_xuanz),1);
for j=1:length(gammaa_xuanz)
    gammaa=gammaa_xuanz(j);
    for i=1:size(tao,1)
        [u_pre(i),s_pre(i)]=eq4(alpha(i), bata, gammaa, tao(i), u0(i), s0(i));
    end
sg=((alpha - bata*u0').*(exp(-bata*tao) - exp(-gammaa*tao)))./(bata - gammaa).^2 + (alpha.*(exp(-gammaa.*tao) - 1))./gammaa.^2 - s0'.*tao.*exp(-gammaa*tao) + (tao.*exp(-gammaa*tao).*(alpha - bata*u0'))./(bata - gammaa) + (alpha.*tao.*exp(-gammaa.*tao))./gammaa;

% f(j)=sum(abs(sg));
f(j)=sum(abs(sg.*(s_obs-s_pre')));
end
[val,indi]=min(f);
gamm=gammaa_xuanz(indi);




end