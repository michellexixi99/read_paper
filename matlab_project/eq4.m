function [ut, st]= eq4(alpha, beta, gamma, tao, u0, s0)
    expu = exp(-beta * tao);
    exps = exp(-gamma * tao);
    ut = u0* expu + alpha / beta * (1 - expu);
    xi = (alpha - u0* beta)/(gamma - beta);
    st = s0* exps + alpha/ gamma * (1 - exps) + xi* (exps - expu);
end
    