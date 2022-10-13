function visualize_model

load('model','theta_f_ests','theta_g_ests','lls','theta_f','theta_g','theta_g_mask');

Dx = length(theta_f);
Dz = size(theta_g,1);
nEpoch = size(theta_f_ests,2)-1;
close all;
figure(1);
cols = {'k' 'r' 'b' 'g' 'c' 'm'};
for i = 1:Dx
    plot(0:nEpoch,theta_f_ests(i,:),[cols{i} '-']); hold on;
    plot([0 nEpoch],[theta_f(i) theta_f(i)],[cols{i} '--']); hold on;
end
ylim([min(theta_f)-0.2, max(theta_f)+0.1]);

figure(2)
for ii = 1:Dx
    subplot(1,Dx,ii)
    idxs = find(theta_g_mask(:,ii)==1);
    for i = 1:length(idxs)
        vec = theta_g_ests(idxs(i),ii,:);
        plot(0:nEpoch,squeeze(vec),[cols{i} '-']); hold on;
        plot([0 nEpoch],[theta_g(idxs(i),ii) theta_g(idxs(i),ii)],[cols{i} '--']); hold on;
    end
end
ylim([min(theta_g(:))-0.1, max(theta_g(:))+0.1]);

figure(100);
subplot(1,2,1)
plot(sum((theta_f_ests-theta_f).^2));
subplot(1,2,2)
mat1 = [];
mat2 = [];
for i = 1:size(theta_g_ests,3)
    m1 = theta_g_ests(:,:,i);
    m2 = theta_g(:,:);
    mat1 = [mat1 m1(:)];
    mat2 = [mat2 m2(:)];
end
plot(sum((mat1-mat2).^2));

figure(200);
plot(lls);
