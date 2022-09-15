function test_wright_fisher_model

N = 100; % size of population
Dz = 3; % dimensionality of Z (# variants)
T = 10;  % # time-points
nSim1 = 20; % # fwd simulations
nSim2 = 10; % # bkw simulations
theta_f = [1  0 -0.5]';  % log relative fitness of variant
theta_h = 0.05 * ones(Dz,1); % mutation rate
theta_z0 = 0.5 * ones(Dz,1); % initial probability of variant
alpha = 0.8; % for fixed alpha
verbose = 0; % verbosity
ep = 0; % for smoothing proposal dist
nEpoch = 10; % # training epochs
    
[Zs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,10,verbose);

P1s = zeros(Dz,T-1);
alphas = zeros(1,T-1);
for i = 1:nSim1
    for t = 1:(T-1)
        P1s(:,t) = P1s(:,t) + squeeze(mean(Zs{i}(t,:,:),2));
        alphas(t) = (length(unique(Pis{i}(t+1,:))) / N);
    end
end
P1s = P1s ./ nSim1;
% alphas = alpha * ones(1,T-1);

% P1s = zeros(1,T);
% for i = 1:nSim1
%     P1s = P1s + mean(Zs{i},2)';
% end
% P1s = P1s ./ nSim1;

prop_Zs = cell(nSim1,nSim2);
prop_Pis = cell(nSim1,nSim2);
prop_log_Qs = zeros(nSim1,nSim2);
prop_log_Ps = zeros(nSim1,nSim2);
for i = 1:nSim1
    
    sd = i*100;
    [Zs_prop, Pis_prop, log_Qs, log_Ps] = wright_fisher_bwd(N,T,nSim2,theta_f,theta_h,theta_z0,sd,verbose,...
        Zs{i}(end,:,:),alphas,P1s,ep);
    
    for j = 1:nSim2
        prop_Zs{i,j} = Zs_prop{j};
        prop_Pis{i,j} = Pis_prop{j};
        prop_log_Qs(i,j) = log_Qs(j);
        prop_log_Ps(i,j) = log_Ps(j);
    end
end

close all;
for ii = 1:Dz
    figure(ii);
    cols = {'k' 'r' 'b' 'g' 'c' 'm'};
    for i = 1:nSim1
        vec = mean(Zs{i}(:,:,ii),2);
        plot(1:T, vec',[cols{mod(i,6)+1} '-'],'linewidth',1.5);
        hold on;
    end
end
for ii = 1:Dz
    figure(100+ii);
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