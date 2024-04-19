function estimate_parameters()
tic

ntasks = feature('numcores'); % number of cpus
poolobj = parpool(ntasks); % start parellel computing

% Define folders
data_dir = './Data/';
output_dir = './Parameter_estimates/';

%% Part 1: estimate regional variables
% Define sampling rate for processing data extraction and storage
fs = 400; % sampling rate data is sent to inference
save_fs = 150; % sampling rate estimates are saved

% Load AAL-atlas brain structure information
load('./Sources/roi_info.mat','selected_roi_index') % load pacellation file
ROIs = unique(cell2mat(selected_roi_index(1,:)));

% Load data
data_path = [data_dir 'example_data.mat']; % can be replaced by your own data (e.g., eeg, meg)
load(data_path,'y')

% Initialise Kalman filter
[xi, ~, A, B, C, Q, R, H, varsigma, v0, N_states, N_syn, N_inputs, N_samples] = forward_simulate_nmm();
Q(10:10, 10:10) = zeros(1, 1); % fix alpha
xi0_empirical = mean(xi(:,N_samples/2:end),2);
P0_empirical = 10*cov(xi(:,N_samples/2:end)');
% P0_empirical(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;
P0_empirical(2*N_syn+1:2*N_syn+1,2*N_syn+1:2*N_syn+1) = eye(N_inputs)*10e-2; % only open "external input"

% Run Kalman filter
for stage = 1:2

if stage == 1
    load([data_dir 'KF_initialization.mat'], 'xi_hat_list')
    % Initialize the parameter vector using estimates obtained from eye-close resting state brain imaging data.
    xi0_array = repmat(xi0_empirical,1,78);
    xi0_array(9,:) = squeeze(mean(xi_hat_list(9,end-2000:end-500,:),[2])); %
    xi0_array(10,:) = squeeze(mean(xi_hat_list(10,end-2000:end-500,:),[2]));
    xi0_array(11,:) = squeeze(mean(xi_hat_list(11,end-2000:end-500,:),[2]));
    xi0_array(12,:) = squeeze(mean(xi_hat_list(12,end-2000:end-500,:),[2]));
    xi0_array(13,:) = squeeze(mean(xi_hat_list(13,end-2000:end-500,:),[2]));
end

if stage == 2
    load([output_dir 'variable_estimates.mat'], 'xi_hat_list')
    % (Optional) run filter second time to make sure it converges
    xi0_array = repmat(xi0_empirical,1,78);
    xi0_array(9,:) = squeeze(max(xi_hat_list(9,1:end-2000,:),[],[2]));
    xi0_array(10,:) = squeeze(max(xi_hat_list(10,1:end-2000,:),[],[2]));
    xi0_array(11,:) = squeeze(max(xi_hat_list(11,1:end-2000,:),[],[2]));
    xi0_array(12,:) = squeeze(max(xi_hat_list(12,1:end-2000,:),[],[2]));
    xi0_array(13,:) = squeeze(max(xi_hat_list(13,1:end-2000,:),[],[2]));
end

clear xi_hat_list

% Inference starts here
parfor (iCh = 1:size(y,1), ntasks)
    measurement = y(iCh,:);      
    xi0 = xi0_array(:,iCh);
    [xi_hat, er] = AKF_quick(measurement, xi0, P0_empirical, A, B, C, Q, R, H, varsigma, v0, N_syn, N_states, N_inputs);
    xi_hat_list(:,:,iCh) = single(resample(xi_hat', save_fs, fs)');
end

% Save regional variables estimates
save([output_dir 'variable_estimates.mat'],'xi_hat_list', '-v7.3');
end

% Display the mean of regional variables estimates with 95% CI
show_estimates_avg_errorband(xi_hat_list);

%% Part 2: estimate inter-regional connectivity
% Define variables before run
window_size_inSeconds = 1; % 1 second time window
save_fs = 150; % sampling rate estimates are saved
fs = 400; % sampling rate data is sent to inference
window_length = window_size_inSeconds*save_fs;

% Define model parameters
v0 = 6;
varsigma = 3.0339;
scale = 50;

%
% load([output_dir 'variable_estimates.mat'],'xi_hat_list');
xi_hat_list(isnan(xi_hat_list))=0; % replace nan with 0
n_channels = size(xi_hat_list,3);
n_timesteps = size(xi_hat_list,2);

xi_hat_list = double(xi_hat_list);
v_pyr = (squeeze(xi_hat_list(1,:,:))' + squeeze(xi_hat_list(7,:,:))' + squeeze(xi_hat_list(9,:,:))')/scale;
v_ip = squeeze(xi_hat_list(1,:,:))'/scale;
z_ip = squeeze(xi_hat_list(2,:,:))';
v_pi = squeeze(xi_hat_list(3,:,:))'/scale; % extract estimates and reshape to [n_channels, n_timesteps]
z_pi = squeeze(xi_hat_list(4,:,:))';
v_pe = squeeze(xi_hat_list(5,:,:))'/scale;
z_pe = squeeze(xi_hat_list(6,:,:))';
v_ep = squeeze(xi_hat_list(7,:,:))'/scale;
z_ep = squeeze(xi_hat_list(8,:,:))';
input = squeeze(xi_hat_list(9,:,:))'/scale;
alpha_ip = fs*squeeze(xi_hat_list(10,:,:))';
alpha_pi = fs*squeeze(xi_hat_list(11,:,:))';
alpha_pe = fs*squeeze(xi_hat_list(12,:,:))';
alpha_ep = fs*squeeze(xi_hat_list(13,:,:))';

% Calculate connection matrix W (mu = W*g(vp))
windows_num = floor(n_timesteps/window_length)-1;

kernel = arrayfun(@h,0:1/save_fs:window_size_inSeconds); % discretized kernel

W_sequence = zeros(n_channels, n_channels, windows_num);
parfor (i = 1:windows_num, ntasks)
    mu = input(:,(i-1)*window_length+1:i*window_length)*scale;
    vp = v_pyr(:,(i-1)*window_length+1:i*window_length);
    gvp = g_function(vp,v0,varsigma);
    
    for j = 1:size(gvp,1) % convolution of the gvp and the kernel
        convolved = conv(gvp(j,:),kernel);
        gvp(j,:) = convolved(1:size(gvp,2));
    end
    
    W_sequence(:,:,i) = analytic_multiregress(gvp', mu');
end
delete(poolobj) % stop parellel computing

% Save inter-regional connectivity estimates
save([output_dir 'connectivity_matrix.mat'], 'W_sequence', '-v7.3')

% Display connectivity estimates
show_connectivity(W_sequence)
toc
end