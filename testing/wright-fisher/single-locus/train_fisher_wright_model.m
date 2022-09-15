function train_fisher_wright_model

N = 100; % size of population FIXED
T = 10;  % # time-points FIXED
nSim1 = 20; % # fwd simulations FIXED training examples
nSim2 = 10; % # bkw simulations FIXED coalescent simulations
theta_f = 0.5;  % log relative fitness of variant INFERRED
theta_h = 0.05; % mutation rate FIXED
theta_z0 = 0.5; % initial probability of variant FIXED
% alpha = 0.8; % for fixed alpha

verbose = 0; % verbosity
ep = 0; % for smoothing proposal dist
nEpoch = 40; % # training epochs

%%% ground truth

sd = 10;
[Zs, Pis, log_Ps] = fisher_wright_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,sd,verbose);

%%% training

theta_f_est = 0;
theta_f_ests = zeros(1,nEpoch+1);
theta_f_ests(1) = theta_f_est;

for cEpoch = 1:nEpoch

    cEpoch
    
    %%% E-step
    % do fwd simulations
    [Zs_fwd, Pis_fwd, log_Ps_fwd] = fisher_wright_fwd(N,T,nSim1,theta_f_est,theta_h,theta_z0,sd*cEpoch*10,verbose);

    prop_Zs = cell(nSim1,nSim2);
    prop_Pis = cell(nSim1,nSim2);
    prop_rs = zeros(nSim1,nSim2);
    for i = 1:nSim1       
        ds = zeros(1,nSim1);
        for j = 1:nSim1
            ds(j) = abs(mean(Zs{i}(end,:)) - mean(Zs_fwd{j}(end,:)));
        end
        fwd_idx = find(ds==min(ds),1);
        alphas = zeros(1,T-1);
        P1s = zeros(1,T-1);
        for t = 1:(T-1)
            P1s(t) = mean(Zs_fwd{fwd_idx}(t,:));
            alphas(t) = (length(unique(Pis_fwd{fwd_idx}(t+1,:))) / N);
%             alphas(t) = alpha;
        end

        sd2 = sd*i*cEpoch*100;
        [Zs_prop, Pis_prop, log_Qs, log_Ps] = fisher_wright_bwd(N,T,nSim2,theta_f_est,theta_h,theta_z0,sd2,verbose,...
            Zs{i}(end,:),alphas,P1s,ep);
        
        for j = 1:nSim2
            prop_Zs{i,j} = Zs_prop{j};
            prop_Pis{i,j} = Pis_prop{j};
        end

        % importance sampling weights
        rs = log_Ps - log_Qs;
        rs = rs - min(rs);        
        rs = exp(rs);
        rs = rs ./ sum(rs);
        prop_rs(i,:) = rs;

    end
    
    %%% M-step
    % update theta_f
    log_ratios = [];
    for i = 1:nSim1
        mat = [];
        for j = 1:nSim2
            vec0 = zeros(T-1,1);
            vec1 = zeros(T-1,1);
            for t = 2:T
                prop1 = mean(prop_Zs{i,j}(t-1,:)==1);
                vec0(t-1) = sum(prop_Zs{i,j}(t-1,prop_Pis{i,j}(t,:))==0)/(1-prop1);
                vec1(t-1) = sum(prop_Zs{i,j}(t-1,prop_Pis{i,j}(t,:))==1)/(prop1);
            end
            vec0 = vec0 + ep;
            vec1 = vec1 + ep;
            mat = [mat (vec1./vec0)];
        end
        log_ratios = [log_ratios ; sum(log(mat) .* repmat(prop_rs(i,:),[T-1 1]),2)];
    end
    theta_f_est = mean(log_ratios(isfinite(log_ratios)));
    theta_f_ests(cEpoch+1) = theta_f_est;
    theta_f_ests(1:cEpoch+1)
    
    if ~isfinite(theta_f_est)
        break;
    end
    
end

close all;
figure(1);
plot(theta_f_ests);
figure(2);
plot((theta_f_ests-theta_f).^2);
