function test_wright_fisher_model

rng(100);
N = 100; % size of population
Dz = 9; % dimensionality of Z (# variants)
Dx = 3; % dimensionality of X (# genes) 
T = 10;  % # time-points
nSim1 = 5; % # fwd simulations
nSim2 = 5; % # bkw simulations
theta_f = [1  0 -0.5]';  % log relative fitness of gene
theta_h = 0.05 * ones(Dz,1); % mutation rate
theta_z0 = 0.5 * ones(Dz,1); % initial probability of variants
theta_g_mask = [repmat([1 0 0],[3 1]) ; repmat([0 1 0],[3 1]); ...
    repmat([0 0 1],[3 1])]; % mask for g-p map
theta_g = randn(Dz,Dx) .* theta_g_mask; % true g-p map
bin_expr_flag = 0; % binary expression flag
% alpha = 0.8; % for fixed alpha
bwd_p_flip = 0.05; % binary flip prob for kl-dist matching
bwd_sigma = 0.1; % gaussian sig value for kl-dist matching
fwd_sigma = 0.1; % gaussian sig value for fwd model
verbose = 0; % verbosity
ep = 0.01; % for smoothing proposal dist
    
sd = 10;
[Zs, Xs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,theta_g,bin_expr_flag,fwd_sigma,sd,verbose);

prop_Zs = cell(nSim1,nSim2);
prop_Xs = cell(nSim1,nSim2);
prop_Pis = cell(nSim1,nSim2);
prop_log_Qs = zeros(nSim1,nSim2);
prop_log_Ps = zeros(nSim1,nSim2);
for i = 1:nSim1
    i
    
    alphas = zeros(1,T-1);
    P1s_Z = zeros(N,Dz,T-1);
    P1s_X = zeros(N,Dx,T-1);
    fwd_idx = ceil(i);
    for t = 1:(T-1)
        P1s_Z(:,:,t) = squeeze(Zs{fwd_idx}(t,:,:));
        P1s_X(:,:,t) = squeeze(Xs{fwd_idx}(t,:,:));
        alphas(t) = length(unique(Pis{fwd_idx}(t+1,:)));
        %             alphas(t) = alpha;
    end
    
    sd2 = i*100;
    [Zs_prop, Xs_prop, Pis_prop, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim2,theta_f,theta_h,theta_z0,...
            theta_g,bin_expr_flag,bwd_p_flip,bwd_sigma,fwd_sigma,...
            sd2,verbose,Zs{i}(end,:,:),Xs{i}(end,:,:),alphas,P1s_Z,P1s_X,ep);
    
    for j = 1:nSim2
        prop_Zs{i,j} = Zs_prop{j};
        prop_Xs{i,j} = Xs_prop{j};
        prop_Pis{i,j} = Pis_prop{j};
        prop_log_Qs(i,j) = log_Qs(j);
        prop_log_Ps(i,j) = log_Ps(j);
    end
end

close all;
figure(1);
for ii = 1:Dz
    subplot(1,Dz,ii);
    cols = {'k' 'r' 'b' 'g' 'c' 'm'};
    for i = 1:nSim1
        vec = mean(Zs{i}(:,:,ii),2);
        plot(1:T, vec',[cols{mod(i,6)+1} '-'],'linewidth',1.5);
        hold on;
    end
end
figure(2);
for ii = 1:Dz
    subplot(1,Dz,ii);
    c = 1;
    for i = 1:nSim1
        for j = 1:nSim2
            vec = mean(prop_Zs{i,j}(:,:,ii),2);
            plot(1:T, vec',[cols{mod(c,6)+1} '-'],'linewidth',1.5);
            hold on;
            c = c + 1;
        end
    end
end
figure(3);
for ii = 1:Dx
    subplot(1,Dx,ii);
    cols = {'k' 'r' 'b' 'g' 'c' 'm'};
    for i = 1:nSim1
        vec = mean(Xs{i}(:,:,ii),2);
        plot(1:T, vec',[cols{mod(i,6)+1} '-'],'linewidth',1.5);
        hold on;
    end
end
figure(4);
for ii = 1:Dx
    subplot(1,Dx,ii);
    c = 1;
    for i = 1:nSim1
        for j = 1:nSim2
            vec = mean(prop_Xs{i,j}(:,:,ii),2);
            plot(1:T, vec',[cols{mod(c,6)+1} '-'],'linewidth',1.5);
            hold on;
            c = c + 1;
        end
    end
end
