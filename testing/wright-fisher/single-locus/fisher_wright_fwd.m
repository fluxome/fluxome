function [Zs, Pis, log_Ps] = fisher_wright_fwd(N,T,nSim,theta_f,theta_h,theta_z0,seed,verbose)

if nargin==0
    N = 100; % size of population
    T = 50;  % # time-points
    nSim = 20; % # simulations
    theta_f = 0.5;  % log relative fitness of variant
    theta_h = 0.05; % mutation rate
    theta_z0 = 0.1; % initial probability of variant
    seed = 100;
end

Zs = cell(1,nSim);
Pis = cell(1,nSim);
log_Ps = zeros(1,nSim);

rng(seed);

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N);
    Pi = zeros(T,N);    
    Z(1,:) = (rand(1,N) < theta_z0);
    log_Ps(cSim) = sum(Z(1,:)*log(theta_z0) + (1-Z(1,:))*log(1-theta_z0));
    
    for t = 2:T
        fs = Z(t-1,:);
        fs = exp(fs .* theta_f);
        fs = fs ./ sum(fs);
        fs_ = cumsum(fs);
        for n = 1:N
            idx = find(rand<=fs_);
            idx = idx(1);
            Pi(t,n) = idx;            
            log_Ps(cSim) = log_Ps(cSim) + log(fs(idx));
            if (rand >= theta_h)
                Z(t,n) = Z(t-1,idx);
                log_Ps(cSim) = log_Ps(cSim) + log(1-theta_h);
            else
                Z(t,n) = 1 - Z(t-1,idx);
                log_Ps(cSim) = log_Ps(cSim) + log(theta_h);
            end
        end
    end
    Zs{cSim} = Z;
    Pis{cSim} = Pi;    
end

if nargin==0
    close all;
    figure(1);
    cols = {'k' 'r' 'b' 'g' 'c' 'm'};
    for i = 1:nSim
        vec = mean(Zs{i},2);
        plot(1:T, vec',[cols{mod(i,6)+1} '-'],'linewidth',1.5);
        hold on;
    end
end
