function tao=solve_tao(alpha,gamma,beta,un,sn,s0,u0)

% alpha=0.83;gamma=9;
% beta=1;
% un=0.377;sn=0.2;s0=0.9;u0=0.9;
    u_inf=alpha/beta;
    s_inf=alpha/gamma;
    betap=beta/(gamma-beta);
    sp=sn-betap*un;
    s_infp=s_inf-betap*u_inf;
    s0p=s0-betap*u0;
    
    if beta>gamma
        tao=-1/gamma*(log((sp-s_infp)/(s0p-s_infp)));
    else
        tao=-1/beta*(log((un-u_inf)/(u0-u_inf)));
    end
tao=max(0.2,tao);

end