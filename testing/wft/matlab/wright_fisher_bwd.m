function [Zs, Xs, Pis, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim,theta_f,theta_h,theta_z0,...
            theta_g,bin_expr_flag,bwd_p_flip,bwd_sigma,fwd_sigma,seed,verbose,Z_T,X_T,alphas,P1s_Z,P1s_X,ep)

Zs = cell(1,nSim);
Xs = cell(1,nSim);
Pis = cell(1,nSim);
log_Qs = zeros(1,nSim);
Dx = length(theta_f);
Dz = size(theta_g,1);

rng(seed);

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N,Dz);
    X = zeros(T,N,Dx);
    Pi = zeros(T,N);
    Z(T,:,:) = Z_T;
    X(T,:,:) = X_T;

    for t = (T-1):-1:1
        P2_Z = squeeze(Z(t+1,:,:));
        P2_X = squeeze(X(t+1,:,:));
        for n = 1:N
            n_sel = ceil(rand*N);
            log_Qs(cSim) = log_Qs(cSim) + log(1/N);
            if (rand >= alphas(t))
                log_Qs(cSim) = log_Qs(cSim) + log(1-alphas(t));
                for i = 1:Dz
                    fl = (rand < bwd_p_flip);
                    Z(t,n,i) = (1-fl)*P1s_Z(n_sel,i,t) + fl*(1-P1s_Z(n_sel,i,t));                    
                    log_Qs(cSim) = log_Qs(cSim) + fl*log(bwd_p_flip) + (1-fl)*log(1-bwd_p_flip);
                end
                if bin_expr_flag
                    for i = 1:Dx
                        fl = (rand < bwd_p_flip);
                        X(t,n,i) = (1-fl)*P1s_X(n_sel,i,t) + fl*(1-P1s_X(n_sel,i,t));
                        log_Qs(cSim) = log_Qs(cSim) + fl*log(bwd_p_flip) + (1-fl)*log(1-bwd_p_flip);
                    end
                else
                    vec = mvnrnd(P1s_X(n_sel,:,t)',bwd_sigma*eye(Dx));
                    X(t,n,:) = vec;
                    [dum val] = mvnpdf_log(vec,P1s_X(n_sel,:,t)',bwd_sigma*eye(Dx));
                    log_Qs(cSim) = log_Qs(cSim) + val;
                end
            else
                log_Qs(cSim) = log_Qs(cSim) + log(alphas(t));
                for i = 1:Dz
                    fl = (rand < bwd_p_flip);
                    Z(t,n,i) = (1-fl)*P2_Z(n_sel,i) + fl*(1-P2_Z(n_sel,i));                    
                    log_Qs(cSim) = log_Qs(cSim) + fl*log(bwd_p_flip) + (1-fl)*log(1-bwd_p_flip);
                end
                if bin_expr_flag
                    for i = 1:Dx
                        fl = (rand < bwd_p_flip);
                        X(t,n,i) = (1-fl)*P2_X(n_sel,i) + fl*(1-P2_X(n_sel,i));
                        log_Qs(cSim) = log_Qs(cSim) + fl*log(bwd_p_flip) + (1-fl)*log(1-bwd_p_flip);
                    end
                else
                    vec = mvnrnd(P2_X(n_sel,:)',bwd_sigma*eye(Dx));
                    X(t,n,:) = vec;
                    [dum val] = mvnpdf_log(vec,P2_X(n_sel,:)',bwd_sigma*eye(Dx));
                    log_Qs(cSim) = log_Qs(cSim) + val;
                end
            end
        end
        fs = zeros(1,N);
        for i = 1:Dx
            fs = fs + X(t,:,i) * theta_f(i);
        end
        fs = exp(fs);        
        fs = fs ./ sum(fs);        
        for n = 1:N
            gs = ones(1,N);
            for i = 1:Dz
                gs(Z(t,:,i)==Z(t+1,n,i)) = gs(Z(t,:,i)==Z(t+1,n,i)) * (1 - theta_h(i));
                gs(Z(t,:,i)~=Z(t+1,n,i)) = gs(Z(t,:,i)~=Z(t+1,n,i)) * theta_h(i);
            end
            preds = squeeze(Z(t,:,:)) * theta_g;
            if bin_expr_flag==1
                preds = 1 ./ (1 + exp(-preds));
                for i = 1:Dx
                    gs = gs .* (X(t,:,i) .* preds(:,i)' + (1-X(t,:,i)) .* (1-preds(:,i)'));
                end
            else
                [ps dum] = mvnpdf_log(squeeze(X(t,:,:)),preds,fwd_sigma*eye(Dx));
                gs = gs .* ps;
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
    Xs{cSim} = X;
    Pis{cSim} = Pi;
end

% eval fwd probs
log_Ps = zeros(1,nSim);
for cSim = 1:nSim
    Z = Zs{cSim};
    X = Xs{cSim};
    Pi = Pis{cSim};
    for i = 1:Dz
        log_Ps(cSim) = log_Ps(cSim) + sum(Z(1,:,i)*log(theta_z0(i)) + (1-Z(1,:,i))*log(1-theta_z0(i)));
    end
    preds = squeeze(Z(1,:,:)) * theta_g;
    if bin_expr_flag==1
        preds = 1 ./ (1 + exp(-preds));
        for i = 1:Dx
            log_Ps(cSim) = log_Ps(cSim) + sum(X(1,:,i).*log(preds(:,i)') + (1-X(1,:,i)).*log(1-preds(:,i)'));
        end
    else
        [dum vec] = mvnpdf_log(squeeze(X(1,:,:)),preds,fwd_sigma*eye(Dx));
        log_Ps(cSim) = log_Ps(cSim) + sum(vec);
    end    
    for t = 2:T
        fs = zeros(1,N);
        for i = 1:Dx
            fs = fs + X(t-1,:,i) * theta_f(i);
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
            preds = squeeze(Z(t,:,:)) * theta_g;
            if bin_expr_flag==1
                preds = 1 ./ (1 + exp(-preds));
                for i = 1:Dx
                    log_Ps(cSim) = log_Ps(cSim) + sum(X(t,:,i).*log(preds(:,i)') + (1-X(t,:,i)).*log(1-preds(:,i)'));
                end
            else
                [dum vec] = mvnpdf_log(squeeze(X(t,:,:)),preds,fwd_sigma*eye(Dx));
                log_Ps(cSim) = log_Ps(cSim) + sum(vec);
            end
        end
    end
end
