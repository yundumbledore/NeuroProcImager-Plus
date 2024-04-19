%% This function is to calculate multiple metrics such as coherence, pdp, rmse, etc
%%
function [output] = compute_metrics(input)
    method = input.method;
    
    if strcmp(method, 'coherence')
        y_prediction = input.y_prediction;
        s_y = input.s_y;
        fs = input.fs;
        figure
        mscohere(y_prediction,s_y,[],[],[],fs)
    elseif strcmp(method, 'converge time')
        xi = input.xi;
        xi_hat = input.xi_hat;
        fs = input.fs;
        [t] = find_last_convergence_time(xi, xi_hat, fs);
        output.convergence_time = t;
    elseif strcmp(method, 'power spectrum density')
        y_prediction = input.y_prediction;
        fs = input.fs;
        figure
        periodogram(y_prediction,rectwin(length(y_predictionp)),length(y_prediction),fs)
    elseif strcmp(method, 'prediction vs signal')
        y_prediction = input.y_prediction;
        s_y = input.s_y;
        figure
        plot(y_prediction)
        figure
        plot(s_y)      
        output = 1;
    elseif strcmp(method, 'parameter estimation vs ground truth')
        xi = input.xi;
        xi_hat = input.xi_hat;
        figure
        plot(xi(9,:))
        hold on
        plot(xi_hat(9,:))
        hold off
        figure
        plot(xi(10,:))
        hold on
        plot(xi_hat(10,:))
        hold off
        figure
        plot(xi(11,:))
        hold on
        plot(xi_hat(11,:))
        hold off
        figure
        plot(xi(12,:))
        hold on
        plot(xi_hat(12,:))
        hold off
        figure
        plot(xi(13,:))
        hold on
        plot(xi_hat(13,:))
        hold off
    elseif strcmp(method, 'prediction rmse')
        y_prediction = input.y_prediction;
        s_y = input.s_y;
%         [y_prediction,~] = envelope(y_prediction,50,'analytic'); % Moving average signal power
%         [s_y,~] = envelope(s_y,134,'analytic'); % Moving average signal power
%         rmse = sqrt(immse(y_prediction, s_y));
        rmse = sqrt(immse(y_prediction, s_y))/(max(s_y) - min(s_y)); % NRMSE oscillation normalised by max min dif
        output.y_rmse = rmse;
    elseif strcmp(method, 'estimation rmse')
        xi = input.xi;
        xi_hat = input.xi_hat;
%         [envelope_up,~] = envelope(xi(1:8,:),50,'analytic'); % Moving average signal power
%         xi(1:8,:) = envelope_up;
%         [envelope_up,~] = envelope(xi_hat(1:8,:),50,'analytic'); % Moving average signal power
%         xi_hat(1:8,:) = envelope_up;
        
%         output.x1_rmse = sqrt(immse(xi(1,:), xi_hat(1,:)));
%         output.x3_rmse = sqrt(immse(xi(3,:), xi_hat(3,:)));
%         output.x5_rmse = sqrt(immse(xi(5,:), xi_hat(5,:)));
%         output.x7_rmse = sqrt(immse(xi(7,:), xi_hat(7,:)));
        
        output.x1_rmse = sqrt(immse(xi(1,:), xi_hat(1,:)))/(max(xi(1,:)) - min(xi(1,:))); % NRMSE of states: oscillation normalised by max min dif
        output.x3_rmse = sqrt(immse(xi(3,:), xi_hat(3,:)))/(max(xi(3,:)) - min(xi(3,:)));
        output.x5_rmse = sqrt(immse(xi(5,:), xi_hat(5,:)))/(max(xi(5,:)) - min(xi(5,:)));
        output.x7_rmse = sqrt(immse(xi(7,:), xi_hat(7,:)))/(max(xi(7,:)) - min(xi(7,:)));
        
        output.x2_rmse = sqrt(immse(xi(2,:), xi_hat(2,:)))/(max(xi(2,:)) - min(xi(2,:))); % NRMSE of 1st order derivatives: oscillation normalised by max min dif
        output.x4_rmse = sqrt(immse(xi(4,:), xi_hat(4,:)))/(max(xi(4,:)) - min(xi(4,:)));
        output.x6_rmse = sqrt(immse(xi(6,:), xi_hat(6,:)))/(max(xi(6,:)) - min(xi(6,:)));
        output.x8_rmse = sqrt(immse(xi(8,:), xi_hat(8,:)))/(max(xi(8,:)) - min(xi(8,:)));
        
%         output.x9_rmse = sqrt(immse(xi(9,:), xi_hat(9,:)));
%         output.x10_rmse = sqrt(immse(xi(10,:), xi_hat(10,:)));
%         output.x11_rmse = sqrt(immse(xi(11,:), xi_hat(11,:)));
%         output.x12_rmse = sqrt(immse(xi(12,:), xi_hat(12,:)));
%         output.x13_rmse = sqrt(immse(xi(13,:), xi_hat(13,:)));
        
        output.x9_rmse = sqrt(immse(xi(9,:), xi_hat(9,:)))/abs(mean(xi(9,:))); % NRMSE constant normalised by itself
        output.x10_rmse = sqrt(immse(xi(10,:), xi_hat(10,:)))/abs(mean(xi(10,:)));
        output.x11_rmse = sqrt(immse(xi(11,:), xi_hat(11,:)))/abs(mean(xi(11,:)));
        output.x12_rmse = sqrt(immse(xi(12,:), xi_hat(12,:)))/abs(mean(xi(12,:)));
        output.x13_rmse = sqrt(immse(xi(13,:), xi_hat(13,:)))/abs(mean(xi(13,:)));
        
        output.states_error = (output.x1_rmse + output.x3_rmse + output.x5_rmse + output.x7_rmse)/4;
        output.parameters_error = (output.x9_rmse + output.x10_rmse + output.x11_rmse...
            + output.x12_rmse + output.x13_rmse)/5;
        output.state_vector_nrmse = (output.x1_rmse + output.x2_rmse + output.x3_rmse + output.x4_rmse + output.x5_rmse +...
            output.x6_rmse + output.x7_rmse + output.x8_rmse + output.x9_rmse + output.x10_rmse + output.x11_rmse +...
            output.x12_rmse + output.x13_rmse)/13;
    end
end

function [t] = find_last_convergence_time(xi, xi_hat, fs)
    Q = [1.5629 0.0156 0.0156 0.0156 0.0156];
    percentage = 0.15;
    xi = xi(9:end,:);
    xi_hat = xi_hat(9:end,:);
    Nparameters = size(xi,1);
    ts = [];
    
    for i = 1:Nparameters % iterate all parameters
        gt = xi(i,1);
        est = xi_hat(i,:);
        sd = Q(i);
        [time] = converge_time(gt, est, percentage, fs, sd);
        ts(i) = time;
    end
    
%     t = ts(1);
    t = max(ts); % find max convergence time as the filter convergence time
    
    if t>28
        t = NaN;
    end
end

function [time] = converge_time(gt, est, percentage, fs, sd)
% To find the convergence time between ground truth and estimation.
%%
    % gt -- ground truth
    % est -- estimation
    % percentage -- percent confidence interval
    % fs -- data sampling rate
    % sd -- standard deviation of ground truth
    
    upper_ci = gt+percentage*abs(gt);
    lower_ci = gt-percentage*abs(gt); 
    
%     upper_ci = gt+percentage*sd;
%     lower_ci = gt-percentage*sd; 
    
    withinCI = (est>lower_ci)&(est<upper_ci);
    reverse_cumsum = cumsum(withinCI,'reverse');
    reference = 1:size(est,2);
    reference = fliplr(reference);
    idx = find(reference == reverse_cumsum,1,'first'); % where filter converges to ground truth
    
    if isempty(idx)
        idx = size(est,2);
    end
    
    time = idx/fs;
end