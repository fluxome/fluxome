function [Zs, Pis, log_Qs, log_Ps] = fisher_wright_bwd(N,T,nSim,theta_f,theta_h,theta_z0,seed,verbose,...
    Z_T,alphas,P1s,ep)

Zs = cell(1,nSim);
Pis = cell(1,nSim);
log_Qs = zeros(1,nSim);

rng(seed);

for t = 1:T-1
    if P1s(t)<ep
        P1s(t) = ep;
    elseif P1s(t)>(1-ep)
        P1s(t) = (1-ep);
    end
end

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N);
    Pi = zeros(T,N);
    Z(T,:) = Z_T;

    for t = (T-1):-1:1
        P2 = mean(Z(t+1,:));
        if P2<ep
            P2 = ep;
        elseif P2>(1-ep)
            P2 = (1-ep);
        end
        for n = 1:N
            if (rand >= alphas(t))
                Z(t,n) = (rand<P1s(t));
                log_Qs(cSim) = log_Qs(cSim) + log(1-alphas(t));
                log_Qs(cSim) = log_Qs(cSim) + Z(t,n)*log(P1s(t)) + (1-Z(t,n))*log(1-P1s(t));
            else
                Z(t,n) = (rand<P2);
                log_Qs(cSim) = log_Qs(cSim) + log(alphas(t));
                log_Qs(cSim) = log_Qs(cSim) + Z(t,n)*log(P2) + (1-Z(t,n))*log(1-P2);                
            end
        end
        fs = Z(t,:);
        fs = exp(fs .* theta_f);
        fs = fs ./ sum(fs);        
        for n = 1:N
            gs = zeros(1,N);
            gs(Z(t,:)==Z(t+1,n)) = (1 - theta_h);
            gs(Z(t,:)~=Z(t+1,n)) = theta_h;
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
    log_Ps(cSim) = sum(Z(1,:)*log(theta_z0) + (1-Z(1,:))*log(1-theta_z0));
    for t = 2:T
        fs = Z(t-1,:);
        fs = exp(fs .* theta_f);
        fs = fs ./ sum(fs);
        for n = 1:N
            idx = Pi(t,n);          
            log_Ps(cSim) = log_Ps(cSim) + log(fs(idx));
            if Z(t,n) == Z(t-1,idx)
                log_Ps(cSim) = log_Ps(cSim) + log(1-theta_h);
            else
                log_Ps(cSim) = log_Ps(cSim) + log(theta_h);
            end
        end
    end
end
