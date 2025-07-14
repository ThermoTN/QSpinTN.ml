function [ tau ] = GetNowtau( Para, beta )
tau = Para.tau_step(end);
for i = 1:length(Para.beta_c)
    if beta < Para.beta_c(i)
        tau = Para.tau_step(i);
        break;
    end
end
end

