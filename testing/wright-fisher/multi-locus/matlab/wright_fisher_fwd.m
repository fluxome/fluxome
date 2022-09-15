function [Zs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim,theta_f,theta_h,theta_z0,seed,verbose)

Zs = cell(1,nSim);
Pis = cell(1,nSim);
log_Ps = zeros(1,nSim);
Dz = length(theta_f);

rng(seed);

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N,Dz);
    Pi = zeros(T,N);
    for i = 1:Dz
        Z(1,:,i) = (rand(1,N) < theta_z0(i));
        log_Ps(cSim) = log_Ps(cSim) + sum(Z(1,:,i)*log(theta_z0(i)) + (1-Z(1,:,i))*log(1-theta_z0(i)));
    end
    
    for t = 2:T
        fs = zeros(1,N);
        for i = 1:Dz
            fs = fs + Z(t-1,:,i) * theta_f(i);
        end
        fs = exp(fs);
        fs = fs ./ sum(fs);
        fs_ = cumsum(fs);
        for n = 1:N
            idx = find(rand<=fs_);
            idx = idx(1);
            Pi(t,n) = idx;            
            log_Ps(cSim) = log_Ps(cSim) + log(fs(idx));
            for i = 1:Dz
                if (rand >= theta_h(i))
                    Z(t,n,i) = Z(t-1,idx,i);
                    log_Ps(cSim) = log_Ps(cSim) + log(1-theta_h(i));
                else
                    Z(t,n,i) = 1 - Z(t-1,idx,i);
                    log_Ps(cSim) = log_Ps(cSim) + log(theta_h(i));
                end
            end
        end
    end
    Zs{cSim} = Z;
    Pis{cSim} = Pi;    
end
