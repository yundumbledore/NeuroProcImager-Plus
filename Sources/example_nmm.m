function [m, er] = example_nmm(y, kf_type, analytic_type, xi0_empirical, P0_empirical, A, B, C, Q, R, H)



%% Options -------------------------------------------------------------- %
% TODO: compare y_prediction, y_real, etc, and figure out what to do with
% that
NStates = 8; % Number of states
NInputs = 1; % Number of external inputs (u)
NParams = 4; % Number of synaptic strength parameters (alpha_ie, alpha_ei, etc...)
NAugmented = NStates + NInputs + NParams; % Total size of augmented state vector

ESTIMATE        = true;         % Run the forward model and estimate (ESTIMATE = true), or just forward (ESTIMATE = false)
PCRB            = 0;            % Compute the PCRB (false = 0, or true > 0) The number here defines the iterations for CRB and MSE
REAL_DATA       = true;         % True to load Seizure activity from neurovista recordings, false to generate data with the forward model
TRUNCATE        = 12000;            % If ~=0, the real data from recordings is truncated from sample 1 to 'TRUNCATE'
SCALE_DATA      = 2;         % Scale Raw data to match dynamic range of the membrane potentials in our model. Multiplies 'y' by the value of SCALE_DATA, try SCALE_DATA = 0.12
INTERPOLATE     = 0;            % Upsample Raw data by interpolating <value> number of samples between each two samples. Doesn't interpolate if INTERPOLATE == {0,1}.

REMOVE_DC       = false;        % Remove DC offset from simulated observed EEG
ADD_NOISE       = true;         % Add noise to the forward model's states
ADD_OBSERVATION_NOISE = true;	% Add noise to the forward model's states

KF_TYPE         = kf_type;  % String: 'unscented', 'extended' (default) or 'none'
ANALYTIC_TYPE   = analytic_type;   % Algorithm to run: 'pip' or 'analytic'. Only makes a difference if the filter (KF_TYPE) is 'extended' or 'none'

ALPHA_KF_LBOUND  = false;       % Zero lower bound (threshold) on alpha in the Kalman Filter (boolean)
ALPHA_KF_UBOUND  = 0;%1e3;      % Upper bound on alpha in the Kalman Filter (integer, if ~=0, the upper bound is ALPHA_KF_UBOUND)
ALPHA_DECAY     = false;        % Exponential decay of alpha-params
FIX_PARAMS      = false;         % Fix input and alpha parameters to initial conditions
RANDOM_ALPHA    = true;         % Chose a random alpha initialization value (true), or the same initialization as the forward model (false)
MONTECARLO      = false;        % Calculatruee term P6 of the covariance matrix (P) by a montecarlo (true), or analytically (false)

relativator = @(x)sqrt(mean(x.^2,2)); % @(x)(max(x')-min(x'))'; % If this is different to @(x)1, it calculates the relative RMSE dividing by whatever this value is.
% ----------------------------------------------------------------------- %

% Initialise random number generator for repeatability
% rng(0);

%% Initialization
% params = set_parameters('alpha', mu); % Set params.u from the input argument 'mu' of set_params
params = set_parameters('alpha');       % Chose params.u from a constant value in set_params

N = 12000; % 148262; % Seizure 1 size: 148262;             	% number of samples
if (TRUNCATE && REAL_DATA), N = TRUNCATE; end % If TRUNCATE ~=0, only take N = TRUNCATE samples of the recording or simulation
% dT = params.dt;         % sampling time step (global)
% dt = 1*dT;            	% integration time step
% nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)
% t = 0:dt:(N-1)*dt;

% model structure
% ~~~~~~~~~~~~~~~~
%           u
%       __  |  ___
%      /  \ | /   \
%    a_ep  \|/    a_ip
%     |     |      |
% ^   E     P      I   ^
% |   |     |      |   |  direction of info
%    a_pe  /|\    a_pi
%     |   / | \    |
%      \_/  |  \__/
%           v
%
% Populations:
%   E - excitatory(pyramidal)
%   I - inhibitory
% Gains:
%   a_ei - connectivity strength from excitatory to inhibitory (alpha_ei)
%   a_ie - connectivity strength from inhibitory to excitatory (alpha_ie)
%

u = params.u;
alpha = [params.alpha_ie; params.alpha_ei; params.alpha_ei; params.alpha_ei]; % Parameters in the augmented state vector. [aip, api, aep, ape]

% Initialise trajectory state
x0 = zeros(NAugmented,1); % initial state
% x0(1:NStates) = mvnrnd(x0(1:NStates),10^1*eye(NStates)); % Random inital state
x0 = params.v0*ones(size(x0));% x0([2 4]) = 0;
x0(NStates+1:end) = [u; alpha];

% Initialize covariance matrix
P0 = 1e2*eye(NAugmented);
% Make P0 different for z-values
P0([2 4 6 8],[2 4 6 8]) = P0([2 4 6 8],[2 4 6 8]) * 50;
% P = zeros(NAugmented, NAugmented, N);
% P(:,:,1) = P0;

% Initialize vector to store firing rates (output of the sigmoid)
% f_i = zeros(1,N); % Firing rate of the inhibitory neurons
% f_e = zeros(1,N); % Firing rate of the excitatory neurons


% Define the model
nmm = nmm_define(x0, P0, params);
% Pass options to model struct
nmm.options.P6_montecarlo   = MONTECARLO;
nmm.options.ALPHA_DECAY     = ALPHA_DECAY;
nmm.options.ALPHA_KF_UBOUND = ALPHA_KF_UBOUND;
nmm.options.ALPHA_KF_LBOUND = ALPHA_KF_LBOUND;
nmm.options.KF_TYPE         = KF_TYPE;
nmm.options.ANALYTIC_TYPE   = ANALYTIC_TYPE;

% Initialize states
% x0 = nmm.x0;
% x = zeros(NAugmented,N);
% x(:,1) = x0;
nmm.A = A;
nmm.B = B;
nmm.C = C;
% Transition model
f = @(x)nmm_run(nmm, x, [], 'transition');
% F = @(x)nmm_run(nmm, x, [], 'jacobian');
% Analytic
f_ = @(x,P)nmm_run(nmm, x, P,  ANALYTIC_TYPE);
% F_ = @(x,P)nmm_run(nmm, x, P,  'jacobian'); % F_ - transition matrix function (a function that takes the current state and returns the Jacobian).

% m0 = mean(xi(:,N_samples/2:end),2);
nmm.x0 = xi0_empirical; % Update initial value in nmm, i.e. nmm.x0
% P_hat_init = 10*cov(xi(:,N_samples/2:end)');
% P_hat_init(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;
nmm.P0 = P0_empirical;

% P0 = 1e2*eye(NAugmented); % P0 will use the same initial value as the
% forward model
% P0(NAugmented - NParams + 1 : end, NAugmented - NParams + 1 : end) = 1e3*eye(NParams);

% Apply EKF filter
try
%     tic
    [m, ~, ~, ~, ~, er] = analytic_kalman_filter_2(y,f_,f,nmm,H,Q,R,'euler');
%     toc_ = toc;
%     disp(['Kalman filter estimation took: ' num2str(toc_) ' seconds']);
catch ME
    if strcmp('Manually stopped', ME.message)
        disp('Kalman filter manually stopped by user.');
        return
    else
        rethrow(ME);
    end
end
end