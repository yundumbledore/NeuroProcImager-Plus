function show_dynamic_stability()
%% Define folders
data_dir = '../Data/';
output_dir = '../Output/';

%% Load data
load([output_dir 'Jacobian_eig.mat'], 'J_eigenvalues_matrix')
load([data_dir 'responsiveness_level.mat'], 'accuracy')
load([data_dir 'gas_level.mat'], 'etconc')

%% Histogram real part of eigenvalues of Jacobi
k = real(J_eigenvalues_matrix);
kk = real(J_eigenvalues_matrix);
kk(k>prctile(k, 99.99, 'all')) = prctile(k, 99.99, 'all');
kk(k<prctile(k, 1, 'all')) = prctile(k, 1, 'all');
clear k J_eigenvalues_matrix

min_val = min(kk,[],'all');
max_val = max(kk,[],'all');
min_val = round(min_val, -1);
max_val = round(max_val, -1);
kk(kk>max_val) = max_val;
kk(kk<min_val) = min_val;

edges = min_val:5:max_val;

for i = 1:size(kk, 2)
    h = histogram(kk(:,i), edges);
    y(:,i) = h.Values;
    close all
end
y = flipud(y);

positive_number = find(edges == max_val) - find(edges == 0);
count = sum(y(1:positive_number,:),1); % find the number of critical mode > 0

%% Adjust count, responsiveness level, gas level for display
startidx = 1;
endidx = numel(count);
countt = count(startidx:endidx);
yy = y(:,startidx:endidx);

accuracyy = accuracy;
countt = movmean(countt, 10);
length = min([numel(countt), numel(accuracyy)]);

%% Plot
fig = figure('Position',[100 100 800 500]);
imagesc(log(yy(:,1:length)))
colormap cool

hold on
plot3(1:length,-1*ones(1,length),countt(1:length),'LineWidth',1)
plot3(1:length,-1*ones(1,length),etconc(1:length),'LineWidth',1)
plot3(1:length,-1*ones(1,length),100*accuracyy(1:length),'LineWidth',1.5)

c = colorbar;
c.Title.String = 'log(count)';
xlabel('Time (seconds)')
ylabel('Real part of eigenvalues')

%% label y ticks
inverse_edges = sort(edges,'descend');
idx40 = find(inverse_edges == 40);
idx20 = find(inverse_edges == 20);
idx0 = find(inverse_edges == 0);
idxneg20 = find(inverse_edges == -20);
idxneg40 = find(inverse_edges == -40);
idxneg50 = find(inverse_edges == -50);

plot3(1:length, idx0*ones(1,length), zeros(1,length),'white','LineWidth',4)
yticks([idx40 idx20 idx0 idxneg20 idxneg40])
yticklabels({'40', '20', '0', '-20', '-40'})
hold off

%% Add legend
legend('Number of positive eigenvalues','Xenon concentration level (%)','Responsiveness level (%)','')

%% set up x, y, z display limit
ylim([-1, idxneg50])

%% set font size
set(gca,'FontSize',12)
end