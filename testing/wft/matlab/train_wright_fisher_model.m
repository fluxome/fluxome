function train_wright_fisher_model

rng(100);

N = 100; % size of population
Dz = 9; % dimensionality of Z (# variants)
Dx = 3; % dimensionality of X (# genes) 
T = 10;  % # time-points
nSim1 = 10; % # fwd simulations
nSim2 = 10; % # bkw simulations
theta_f = [1  0 -0.5]';  % log relative fitness of gene
theta_h = 0.05 * ones(Dz,1); % mutation rate
theta_z0 = 0.5 * ones(Dz,1); % initial probability of variants
theta_g_mask = [repmat([1 0 0],[3 1]) ; repmat([0 1 0],[3 1]); ...
    repmat([0 0 1],[3 1])]; % mask for g-p map
theta_g = randn(Dz,Dx) .* theta_g_mask; % true g-p map
bin_expr_flag = 1; % binary expression flag
% alpha = 0.8; % for fixed alpha
bwd_p_flip = 0.05; % binary flip prob for kl-dist matching
bwd_sigma = 0.5; % gaussian sig value for kl-dist matching
fwd_sigma = 0.2; % gaussian sig value for fwd model
verbose = 0; % verbosity
ep = 0.01; % for smoothing proposal dist
nEpoch = 10; % # training epochs

%%% ground truth

sd = 10;
[Zs, Xs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,theta_g,bin_expr_flag,fwd_sigma,sd,verbose);

%%% training

theta_f_est = zeros(Dx,1);
theta_f_ests = zeros(Dx,nEpoch+1);
theta_f_ests(:,1) = theta_f_est;
theta_g_est = zeros(Dz,Dx,1);
theta_g_ests = zeros(Dz,Dx,nEpoch+1);
theta_g_ests(:,:,1) = theta_g_est;
lls = [];

for cEpoch = 1:nEpoch

    cEpoch
    
    %%% E-step
    % do fwd simulations
    [Zs_fwd, Xs_fwd, Pis_fwd, log_Ps_fwd] = wright_fisher_fwd(N,T,nSim1,theta_f_est,theta_h,theta_z0,...
        theta_g_est,bin_expr_flag,fwd_sigma,sd*cEpoch*10,verbose);

    prop_Zs = cell(nSim1,nSim2);
    prop_Xs = cell(nSim1,nSim2);
    prop_Pis = cell(nSim1,nSim2);
    prop_rs = zeros(nSim1,nSim2);
    for i = 1:nSim1     
        [cEpoch i]
        ds = zeros(1,nSim1);
        for j = 1:nSim1
            
            %%% add kl test
            log_PZ_mat = zeros(N,N);
            for n = 1:N
                log_PZ_mat1 = (squeeze(Zs{i}(end,:,:)) - repmat(squeeze(Zs_fwd{j}(end,n,:))',[N 1])).^2;
                log_PZ_mat1 = log(bwd_p_flip) * log_PZ_mat1 + log(1-bwd_p_flip) * (1-log_PZ_mat1);
                log_PZ_mat(:,n) = sum(log_PZ_mat1,2);
            end
            if bin_expr_flag == 1
                log_PX_mat = zeros(N,N);
                for n = 1:N
                    log_PX_mat1 = (squeeze(Xs{i}(end,:,:)) - repmat(squeeze(Xs_fwd{j}(end,n,:))',[N 1])).^2;
                    log_PX_mat1 = log(bwd_p_flip) * log_PX_mat1 + log(1-bwd_p_flip) * (1-log_PX_mat1);
                    log_PX_mat(:,n) = sum(log_PX_mat1,2);
                end
            else
                log_PX_mat = zeros(N,N);
                for n = 1:N
                    [dum vec] = mvnpdf_log(squeeze(squeeze(Xs{i}(end,:,:)), Xs_fwd{j}(end,n,:)), bwd_sigma*eye(Dx));
                    log_PX_mat(:,n) = vec;
                end                    
            end
            ds(j) = sum(logsumexp((log(1/N) + log_PZ_mat + log_PX_mat),2));
        end
        fwd_idx = find(ds==min(ds),1);
        alphas = zeros(1,T-1);
        P1s_Z = zeros(N,Dz,T-1);
        P1s_X = zeros(N,Dx,T-1);
        for t = 1:(T-1)
            P1s_Z(:,:,t) = squeeze(Zs_fwd{fwd_idx}(t,:,:));
            P1s_X(:,:,t) = squeeze(Xs_fwd{fwd_idx}(t,:,:));
            alphas(t) = (length(unique(Pis_fwd{fwd_idx}(t+1,:))) / N);
%             alphas(t) = alpha;
        end

        sd2 = sd*i*cEpoch*100;
        [Zs_prop, Xs_prop, Pis_prop, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim2,theta_f_est,theta_h,theta_z0,...
            theta_g_est,bin_expr_flag,bwd_p_flip,bwd_sigma,fwd_sigma,...
            sd2,verbose,Zs{i}(end,:,:),Xs{i}(end,:,:),alphas,P1s_Z,P1s_X,ep);
        
        for j = 1:nSim2
            prop_Zs{i,j} = Zs_prop{j};
            prop_Xs{i,j} = Xs_prop{j};            
            prop_Pis{i,j} = Pis_prop{j};
        end
        rs = log_Ps - log_Qs;
        rs = rs - max(rs);
        rs = exp(rs);
        rs = rs ./ sum(rs);
        prop_rs(i,:) = rs;
    end
    
    %%% M-step
    % update theta_f
    fun = @(x)-loglike(prop_Zs,prop_Xs,prop_Pis,prop_rs,theta_g_mask,bin_expr_flag,fwd_sigma,x);
    opts = optimoptions('fminunc','MaxIterations',2,'Display','iter');
    [theta_est_new,fval] = fminunc(fun,[theta_f_est ; theta_g_est(theta_g_mask(:)~=0)],opts);    
    theta_f_est = theta_est_new(1:Dx);
    theta_f_ests(:,cEpoch+1) = theta_f_est;
    theta_g_est(theta_g_mask(:)~=0) = theta_est_new(Dx+1:end);
    theta_g_ests(:,:,cEpoch+1) = theta_g_est;    
    lls = [lls -fval];
    lls
    
    if ~isfinite(theta_f_est)
        break;
    end

    save('model','theta_f_ests','theta_g_ests','lls','theta_f','theta_g','theta_g_mask','bin_expr_flag');
    
end

close all;
figure(1);
cols = {'k' 'r' 'b' 'g' 'c' 'm'};
for i = 1:Dx
    plot(0:nEpoch,theta_f_ests(i,:),[cols{i} '-']); hold on;
    plot([0 nEpoch],[theta_f(i) theta_f(i)],[cols{i} '--']); hold on;
end
ylim([min(theta_f)-0.1, max(theta_f)+0.1]);
figure(100);
plot(sum((theta_f_ests-theta_f).^2));
figure(200);
plot(lls);

%%%%% Opt function

function val = loglike(Zs,Xs,Pis,rs,theta_g_mask,bin_expr_flag,fwd_sigma,theta)

T = size(Zs{1,1},1);
N = size(Zs{1,1},2);
Dz = size(Zs{1,1},3);
Dx = size(Xs{1,1},3);
nSim1 = size(Zs,1);
nSim2 = size(Zs,2);

theta_f = theta(1:Dx);
theta_g = zeros(Dz,Dx);
theta_g(theta_g_mask(:)~=0) = theta(Dx+1:end);

val = 0;

for i = 1:nSim1
    for j = 1:nSim2
        for t = 1:T-1
            for n = 1:N
                M = sum(Pis{i,j}(t+1,:)==n);
                val = val + M * sum(theta_f.*squeeze(Xs{i,j}(t,n,:))) * rs(i,j);
            end
        end
        for t = 1:T-1
            val1 = 0;
            for n = 1:N
                val1 = val1 + exp(sum(theta_f.*squeeze(Xs{i,j}(t,n,:))));
            end
            val = val - N * log(val1) * rs(i,j);
        end
        for t = 1:T
            preds = squeeze(Zs{i}(t,:,:)) * (theta_g .* theta_g_mask);
            if bin_expr_flag==1
                preds = 1 ./ (1 + exp(-preds));
                mat = squeeze(Xs{i}(t,:,:)) .* log(preds) + (1-squeeze(Xs{i}(t,:,:))) .* log(1-preds);
                val = val + sum(mat(:)) * rs(i,j);
            else
                [dum vec] = mvnpdf_log(squeeze(Xs{i}(t,:,:)),preds,fwd_sigma*eye(Dx));
                val = val + sum(vec) * rs(i,j);
            end
        end   
    end
end

