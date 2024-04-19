function show_estimates_avg_errorband(xi_hat_list)
    ACPT_dir = '../Data/responsiveness_level.mat'; 
    GAS_dir = '../Data/gas_level.mat'; 

    load(ACPT_dir, 'accuracy') % load responsiveness level data
    load(GAS_dir, 'etconc','t') % load xenon concentration level data
    accuracy = accuracy*100;

    xi_hat_list = xi_hat_list(:,150*10:end,:);
    
    % subplot shows all variable estimates
    figure('Position',[100 100 1200 700]);
    subplot(2,3,1);
    yyaxis left
    y = nanmean(xi_hat_list(9,:,:),[3]);
    stderr = nanstd(squeeze(xi_hat_list(9,:,:)),0,2)'/sqrt(78);
    x=plot_mean_errorband(y,stderr);
    ylabel('External input')
    yyaxis right
    plot(1:numel(accuracy),accuracy)
    hold on
    plot(t,etconc)
    hold off
    ylim([-1 101])
    ylabel('Responsiveness level & Xenon concentration level')
    
    subplot(2,3,2);
    yyaxis left
    y = nanmean(xi_hat_list(10,:,:),[3]);
    stderr = nanstd(squeeze(xi_hat_list(10,:,:)),0,2)'/sqrt(78);
    x=plot_mean_errorband(y,stderr);
    ylabel('aIP')
    yyaxis right
    plot(1:numel(accuracy),accuracy)
    hold on
    plot(t,etconc)
    hold off
    ylim([-1 101])
    ylabel('Responsiveness level & Xenon concentration level')
    
    subplot(2,3,3);
    yyaxis left
    y = nanmean(xi_hat_list(11,:,:),[3]);
    stderr = nanstd(squeeze(xi_hat_list(11,:,:)),0,2)'/sqrt(78);
    x=plot_mean_errorband(y,stderr);
    ylabel('aPI')
    yyaxis right
    plot(1:numel(accuracy),accuracy)
    hold on
    plot(t,etconc)
    hold off
    ylim([-1 101])
    ylabel('Responsiveness level & Xenon concentration level')
    
    subplot(2,3,4);
    yyaxis left
    y = nanmean(xi_hat_list(12,:,:),[3]);
    stderr = nanstd(squeeze(xi_hat_list(12,:,:)),0,2)'/sqrt(78);
    x=plot_mean_errorband(y,stderr);
    ylabel('aPE')
    yyaxis right
    plot(1:numel(accuracy),accuracy)
    hold on
    plot(t,etconc)
    hold off
    ylim([-1 101])
    ylabel('Responsiveness level & Xenon concentration level')
    
    subplot(2,3,5);
    yyaxis left
    y = nanmean(xi_hat_list(13,:,:),[3]);
    stderr = nanstd(squeeze(xi_hat_list(13,:,:)),0,2)'/sqrt(78);
    x=plot_mean_errorband(y,stderr);
    ylabel('aEP')
    yyaxis right
    plot(1:numel(accuracy),accuracy)
    hold on
    plot(t,etconc)
    hold off
    ylim([-1 101])
    ylabel('Responsiveness level & Xenon concentration level')

    % add legend
    Lgnd = legend('95% errorband of mean of estimates', 'Mean of estimates across all regions',...
        'Responsiveness level (%)', 'Xenon concentration level (%)');
    Lgnd.Position(1) = 0.7;
    Lgnd.Position(2) = 0.4;
    Lgnd.FontSize = 13;

    % add title
    sgtitle('Regional Variables Estimates Mean with 95% Errorband across All Regions')
end

function x = plot_mean_errorband(y,stderr)
    x = (1:numel(y))/150;
    curve1 = y + 1.96*stderr;
    curve2 = y - 1.96*stderr;
    
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    fill(x2, inBetween, 'g','LineStyle','none');
    hold on;
    plot(x, y);
    xlabel('Time (Seconds)')
    hold off
end