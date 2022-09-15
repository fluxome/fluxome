function test_fisher_wright_model

N = 100; % size of population
T = 50;  % # time-points
nSim1 = 20; % # fwd simulations
nSim2 = 5; % # bkw simulations
theta_f = 0.5;  % log relative fitness of variant
theta_h = 0.05; % mutation rate
theta_z0 = 0.1; % initial probability of variant
alpha = 0.85*ones(1,T); % mixing param for bkw simulation
verbose = 0; % verbosity
ep = 0; % for smoothing proposal dist
    
[Zs, Pis, log_Ps] = fisher_wright_fwd(N,T,nSim1,theta_f,theta_h,theta_z0,10,verbose);

P1s = zeros(1,T);
for i = 1:nSim1
    P1s = P1s + mean(Zs{i},2)';
end
P1s = P1s ./ nSim1;

prop_Zs = cell(nSim1,nSim2);
prop_Pis = cell(nSim1,nSim2);
prop_log_Qs = zeros(nSim1,nSim2);
prop_log_Ps = zeros(nSim1,nSim2);
for i = 1:nSim1
    
    sd = i*100;
    [Zs_prop, Pis_prop, log_Qs, log_Ps] = fisher_wright_bwd(N,T,nSim2,theta_f,theta_h,theta_z0,sd,verbose,...
        Zs{i}(end,:),alpha,P1s,ep);
    
    for j = 1:nSim2
        prop_Zs{i,j} = Zs_prop{j};
        prop_Pis{i,j} = Pis_prop{j};
        prop_log_Qs(i,j) = log_Qs(j);
        prop_log_Ps(i,j) = log_Ps(j);
    end
end

close all;
figure(1);
cols = {'k' 'r' 'b' 'g' 'c' 'm'};
for i = 1:nSim1
    vec = mean(Zs{i},2);
    plot(1:T, vec',[cols{mod(i,6)+1} '-'],'linewidth',1.5);
    hold on;
end
figure(2);
c = 1;
for i = 1:nSim1
    for j = 1:nSim2
        vec = mean(prop_Zs{i,j},2);
        plot(1:T, vec',[cols{mod(c,6)+1} '-'],'linewidth',1.5);
        hold on;
        c = c + 1;
    end
end
