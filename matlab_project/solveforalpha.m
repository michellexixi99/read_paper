function alphaa=solveforalpha(u_obs,u_pre,s_obs,s_pre,u0,gamma,t, bata,tao, s0)
    t0=t'-tao;
alpha_xuanz=0.01:0.01:1;
for j=1:length(alpha_xuanz)
    alpha=alpha_xuanz(j);
    for i=1:size(u_obs,1)
        [u_pre(i),s_pre(i)]=eq4(alpha, bata, gamma, tao(i), u0(i), s0(i));
    end
    
    ua=-(exp(-bata*(tao)) - 1)/bata;
    sa=(exp(-bata*(tao)) - exp(-gamma*(tao)))/(bata - gamma) - (exp(-gamma*(tao)) - 1)/gamma;
    f(j)=abs(sum(ua*sum(u_obs-u_pre)+sa*sum(s_obs-s_pre)));
end

[val,indi]=min(f);
alphaa=alpha_xuanz(indi);


end