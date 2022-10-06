function [Zs, Pis, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim,theta_f,theta_h,theta_z0,seed,verbose,...
    Z_T,alphas,P1s,ep)

Zs = cell(1,nSim);
Pis = cell(1,nSim);
log_Qs = zeros(1,nSim);
Dz = length(theta_f);

rng(seed);

for t = 1:T-1
    for i = 1:Dz
        if P1s(i,t)<ep
            P1s(i,t) = ep;
        elseif P1s(i,t)>(1-ep)
            P1s(i,t) = (1-ep);
        end
    end
end

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N,Dz);
    Pi = zeros(T,N);
    Z(T,:,:) = Z_T;

    for t = (T-1):-1:1
        P2 = squeeze(mean(Z(t+1,:,:),2));
        for i = 1:Dz
            if P2(i)<ep
                P2(i) = ep;
            elseif P2(i)>(1-ep)
                P2(i) = (1-ep);
            end
        end
        for n = 1:N
            if (rand >= alphas(t))
                log_Qs(cSim) = log_Qs(cSim) + log(1-alphas(t));     % altered
                for i = 1:Dz
                    Z(t,n,i) = (rand<P1s(i,t));
                    log_Qs(cSim) = log_Qs(cSim) + Z(t,n,i)*log(P1s(i,t)) + (1-Z(t,n,i))*log(1-P1s(i,t));
                end
            else
                log_Qs(cSim) = log_Qs(cSim) + log(alphas(t));       % altered
                for i = 1:Dz
                    Z(t,n,i) = (rand<P2(i));                 
                    log_Qs(cSim) = log_Qs(cSim) + Z(t,n,i)*log(P2(i)) + (1-Z(t,n,i))*log(1-P2(i));
                end
            end
        end
        fs = zeros(1,N);
        for i = 1:Dz
            fs = fs + Z(t,:,i) * theta_f(i);
        end
        fs = exp(fs);        
        fs = fs ./ sum(fs);        
        for n = 1:N
            gs = ones(1,N);
            for i = 1:Dz
                gs(Z(t,:,i)==Z(t+1,n,i)) = gs(Z(t,:,i)==Z(t+1,n,i)) * (1 - theta_h(i));
                gs(Z(t,:,i)~=Z(t+1,n,i)) = gs(Z(t,:,i)~=Z(t+1,n,i)) * theta_h(i);
            end
            fgs = fs.*gs;
            fgs = fgs ./ sum(fgs);
            fgs_ = cumsum(fgs);    
            idx = find(rand<=fgs_);
            idx = idx(1);
            Pi(t+1,n) = idx;
            log_Qs(cSim) = log_Qs(cSim) + log(fgs(idx));
        end
    end
    Zs{cSim} = Z;
    Pis{cSim} = Pi;
end

% eval fwd probs
log_Ps = zeros(1,nSim);
for cSim = 1:nSim
    Z = Zs{cSim};
    Pi = Pis{cSim};
    for i = 1:Dz
        log_Ps(cSim) = log_Ps(cSim) + sum(Z(1,:,i)*log(theta_z0(i)) + (1-Z(1,:,i))*log(1-theta_z0(i)));
    end
    for t = 2:T
        fs = zeros(1,N);
        for i = 1:Dz
            fs = fs + Z(t-1,:,i) * theta_f(i);
        end
        fs = exp(fs);        
        fs = fs ./ sum(fs);
        for n = 1:N
            idx = Pi(t,n);          
            log_Ps(cSim) = log_Ps(cSim) + log(fs(idx));
            for i = 1:Dz
                if Z(t,n,i) == Z(t-1,idx,i)
                    log_Ps(cSim) = log_Ps(cSim) + log(1-theta_h(i));
                else
                    log_Ps(cSim) = log_Ps(cSim) + log(theta_h(i));
                end
            end
        end
    end
end
