% [620 -950 260 1100 900] default filter initial parameters values
% from estimation: aEP 112; aIP -70; aPE 720; aPI 757; input 383

tic
final_results = zeros(1000, 32);
for rng_i = 1:30 % introduce random number seeds
%%
ntasks = 36;
fs = 400;
it = 20;
input = 640:it:820;
aIP = -1030:it:-850;
aPI = 210; 
aPE = 1050;
aEP = 820:it:1000; 

saved_file_name = './Results/parameter_space_explore_30iter.csv';
%%
[a,b,c,d,e] = ndgrid(input,aIP,aPI,aPE,aEP);
Z = [a(:),b(:),c(:),d(:),e(:)]; 

input = 300;
input_offset = [];

% parameters = Z(1,:);
% parameters = [600 -1000 200 1000 950];
Niter = size(Z, 1);
results = zeros(size(Z,1), 26);
peakfreqs = zeros(size(Z,1),1);

poolobj = parpool(ntasks);
parfor (i = 1:Niter, ntasks)
%     lineLength = fprintf('Running...%i/%i',i,Niter);
    %%
    parameters = Z(i,:);
    parameters = transpose(parameters).*ones(5,12000);
    [s_y,A, B, C, phi, xi] = generate_simulated_data(input, input_offset, parameters, fs, rng_i);
%     s_y = s_y - mean(s_y);% fs = 400;
    
    [fmax] = find_peak_frequency(s_y, fs) % Check peak frequency
    peakfreqs(i,1) = fmax;
    %%
    [outarg] = diy(s_y, xi, fs); % Run akf and ukf and compare performance
    results(i,:) = outarg;
    %%
%     fprintf(repmat('\b',1,lineLength))
end
delete(poolobj)

estimation_performance = [Z peakfreqs results]; % Z column 1-5 parameters ground truth

final_results = final_results + estimation_performance;
% save('./Results/Z_3parameters_within200_new.mat','Z')
end
% peakfreqs = [Z peakfreqs];
% save('./Results/parameters_space_explore.mat','estimation_performance','-v7.3')
estimation_performance = final_results/30;

T = array2table(estimation_performance);
T.Properties.VariableNames(1:32) = {'input','aIP','aPI','aPE','aEP',...
    'peak_frequency',...
    'akf_prediction_rmse','akf_vip_rmse','akf_vpi_rmse','akf_vpe_rmse','akf_vep_rmse','akf_input_rmse','akf_aIP_rmse','akf_aPI_rmse','akf_aPE_rmse','akf_aEP_rmse','akf_states_rmse','akf_parameters_rmse','akf_convergence_time',...
    'ukf_prediction_rmse','ukf_vip_rmse','ukf_vpi_rmse','ukf_vpe_rmse','ukf_vep_rmse','ukf_input_rmse','ukf_aIP_rmse','ukf_aPI_rmse','ukf_aPE_rmse','ukf_aEP_rmse','ukf_states_rmse','ukf_parameters_rmse','ukf_convergence_time'};
writetable(T, saved_file_name)

%%
toc
function [fmax] = find_peak_frequency(s_y, fs)
    L = 2000;
    noverlap = 1000;
    nfft = 2000;
    
    [pxx,f] = pwelch(s_y, L, noverlap, nfft, fs);
    
    [~,B] = max(pxx);
    fmax = f(B);
end

function [outarg] = diy(s_y, xi, fs)
    outarg = zeros(1,26);
    kf_type = 'extended'; % 'unscented', 'extended'
    analytic_type = 'pip'; % 'pip'
    [y_prediction, xi_hat, er] = example_nmm(s_y, kf_type, analytic_type, fs);
    if er
        outarg(1,:) = nan(1,26);
        return
    end
    
    inputarg.y_prediction = y_prediction;%(end - fs*5 + 1:end); % only take the last 5 seconds
    inputarg.s_y = s_y;%(end - fs*5 + 1:end);
    inputarg.method = 'prediction rmse';
    [output] = compute_metrics(inputarg);
    outarg(1, 1) = output.y_rmse;
    clear inputarg output
    
    inputarg.xi = xi;%(:,end - fs*5 + 1:end);
    inputarg.xi_hat = xi_hat;%(:,end - fs*5 + 1:end);
    inputarg.method = 'estimation rmse';
    [output] = compute_metrics(inputarg);
    outarg(1, 2) = output.x1_rmse;
    outarg(1, 3) = output.x3_rmse;
    outarg(1, 4) = output.x5_rmse;
    outarg(1, 5) = output.x7_rmse;
    outarg(1, 6) = output.x9_rmse;
    outarg(1, 7) = output.x10_rmse;
    outarg(1, 8) = output.x11_rmse;
    outarg(1, 9) = output.x12_rmse;
    outarg(1, 10) = output.x13_rmse;
    outarg(1, 11) = output.states_error;
    outarg(1, 12) = output.parameters_error;
    clear inputarg output
    
    inputarg.xi = xi;
    inputarg.xi_hat = xi_hat;
    inputarg.fs = fs;
    inputarg.method = 'converge time';
    [output] = compute_metrics(inputarg);
    outarg(1, 13) = output.convergence_time;
    
    %%
    kf_type = 'unscented'; % 'unscented', 'extended'
    analytic_type = 'pip'; % 'pip'
    [y_prediction, xi_hat, er] = example_nmm(s_y, kf_type, analytic_type, fs);
     if er
        outarg(1,:) = nan(1,26);
        return
    end
    
    inputarg.y_prediction = y_prediction;%(end - fs*5 + 1:end);
    inputarg.s_y = s_y;%(end - fs*5 + 1:end);
    inputarg.method = 'prediction rmse';
    [output] = compute_metrics(inputarg);
    outarg(1, 14) = output.y_rmse;
    clear inputarg output
    
    inputarg.xi = xi;%(:,end - fs*5 + 1:end);
    inputarg.xi_hat = xi_hat;%(:,end - fs*5 + 1:end);
    inputarg.method = 'estimation rmse';
    [output] = compute_metrics(inputarg);
    outarg(1, 15) = output.x1_rmse;
    outarg(1, 16) = output.x3_rmse;
    outarg(1, 17) = output.x5_rmse;
    outarg(1, 18) = output.x7_rmse;
    outarg(1, 19) = output.x9_rmse;
    outarg(1, 20) = output.x10_rmse;
    outarg(1, 21) = output.x11_rmse;
    outarg(1, 22) = output.x12_rmse;
    outarg(1, 23) = output.x13_rmse;
    outarg(1, 24) = output.states_error;
    outarg(1, 25) = output.parameters_error;
    clear inputarg output
    
    inputarg.xi = xi;
    inputarg.xi_hat = xi_hat;
    inputarg.fs = fs;
    inputarg.method = 'converge time';
    [output] = compute_metrics(inputarg);
    outarg(1, 26) = output.convergence_time;
end





