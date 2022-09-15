function visualize_model

load('model','theta_f_ests','lls','theta_f');

Dz = length(theta_f);
nEpoch = size(theta_f_ests,2)-1;
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
