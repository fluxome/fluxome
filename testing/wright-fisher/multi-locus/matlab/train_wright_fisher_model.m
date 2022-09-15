function train_wright_fisher_model

rng(100);

N = 100; % size of population
Dz = 3; % dimensionality of Z (# variants)
T = 10;  % # time-points
nSim1 = 20; % # fwd simulations
nSim2 = 10; % # bkw simulations
theta_f = [1  0 -0.5]';  % log relative fitness of variant
theta_h = 0.05 * ones(Dz,1); % mutation rate
theta_z0 = 0.5 * ones(Dz,1); % initial probability of variant
% alpha = 0.8; % for fixed alpha
verbose = 0; % verbosity
ep = 0.01; % for smoothing proposal dist
nEpoch = 10; % # training epochs

%%% ground truth

sd = 10;
[Zs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,sd,verbose);

%%% training

theta_f_est = zeros(Dz,1);
theta_f_ests = zeros(Dz,nEpoch+1);
theta_f_ests(:,1) = theta_f_est;
lls = [];

for cEpoch = 1:nEpoch

    cEpoch
    
    %%% E-step
    % do fwd simulations
    [Zs_fwd, Pis_fwd, log_Ps_fwd] = wright_fisher_fwd(N,T,nSim1,theta_f_est,theta_h,theta_z0,sd*cEpoch*10,verbose);

    prop_Zs = cell(nSim1,nSim2);
    prop_Pis = cell(nSim1,nSim2);
    prop_rs = zeros(nSim1,nSim2);
    for i = 1:nSim1       
        ds = zeros(1,nSim1);
        for j = 1:nSim1
            ds(j) = sum((mean(Zs{i}(end,:,:),2) - mean(Zs_fwd{j}(end,:,:),2)).^2);
        end
        fwd_idx = find(ds==min(ds),1);
        alphas = zeros(1,T-1);
        P1s = zeros(Dz,T-1);
        for t = 1:(T-1)
            P1s(:,t) = mean(Zs_fwd{fwd_idx}(t,:,:),2);
            alphas(t) = (length(unique(Pis_fwd{fwd_idx}(t+1,:))) / N);
%             alphas(t) = alpha;
        end

        sd2 = sd*i*cEpoch*100;
        [Zs_prop, Pis_prop, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim2,theta_f_est,theta_h,theta_z0,sd2,verbose,...
            Zs{i}(end,:,:),alphas,P1s,ep);
        
        for j = 1:nSim2
            prop_Zs{i,j} = Zs_prop{j};
            prop_Pis{i,j} = Pis_prop{j};
        end
        rs = log_Ps - log_Qs;
        rs = rs - min(rs);
        if i==1
%             rs
        end
        rs = exp(rs);
        rs = rs ./ sum(rs);
        prop_rs(i,:) = rs;
    end
    
    %%% M-step
    % update theta_f
    theta_f_est
    fun = @(x)-loglike(prop_Zs,prop_Pis,prop_rs,x);
    opts = optimoptions('fminunc','MaxIterations',2,'Display','iter');
    [theta_f_est_new,fval] = fminunc(fun,theta_f_est,opts);
    theta_f_est = theta_f_est_new;
    theta_f_ests(:,cEpoch+1) = theta_f_est;
    lls = [lls -fval];
    lls
    
    if ~isfinite(theta_f_est)
        break;
    end

    save('model','theta_f_ests','lls','theta_f');
    
end

close all;
figure(1);
cols = {'k' 'r' 'b' 'g' 'c' 'm'};
for i = 1:Dz
    plot(0:nEpoch,theta_f_ests(i,:),[cols{i} '-']); hold on;
    plot([0 nEpoch],[theta_f(i) theta_f(i)],[cols{i} '--']); hold on;
end
ylim([min(theta_f)-0.1, max(theta_f)+0.1]);
figure(100);
plot(sum((theta_f_ests-theta_f).^2));
figure(200);
plot(lls);

%%%%% Opt function

function val = loglike(Zs,Pis,rs,theta_f)

T = size(Zs{1,1},1);
N = size(Zs{1,1},2);
nSim1 = size(Zs,1);
nSim2 = size(Zs,2);

val = 0;

for i = 1:nSim1
    for j = 1:nSim2
        for t = 1:T-1
            for n = 1:N
                M = sum(Pis{i,j}(t+1,:)==n);
                val = val + M * sum(theta_f.*squeeze(Zs{i,j}(t,n,:))) * rs(i,j);
            end
        end
        for t = 1:T-1
            val1 = 0;
            for n = 1:N
                val1 = val1 + exp(sum(theta_f.*squeeze(Zs{i,j}(t,n,:))));
            end
            val = val - N * log(val1) * rs(i,j);
        end
    end
end

