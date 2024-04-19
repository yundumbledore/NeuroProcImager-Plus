% This code is to do kalman filter on local machine
function [xi_hat,v_pyr_hat,xii]=kalman_estimation_2(s_y,fs,scaler,rr)
    % generate some model data
%     disp(['Processing channel ' num2str(ich) '...'])
    
    input = 300;
    input_offset = [];

    [A,B,C,N_states,N_syn,N_inputs,N_samples,xi,v0,varsigma,Q,R,H] = set_params(input,input_offset,fs,scaler,rr);
    xii = xi;
%     save('~/Desktop/initials.mat','A','B','C','N_states','N_syn','N_inputs','N_samples','xi','v0','varsigma','Q','R','H')

    xi_hat_init = mean(xi(:,N_samples/2:end),2);                                % to ignore inital transient take the mean of the second half of the test data
    P_hat_init = 10*cov(xi(:,N_samples/2:end)');
    P_hat_init(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;               % open up the error in connectivity strength.
    
    xi_hat = zeros(N_states,numel(s_y)+1);                 % initialize for speed
    P_hat = zeros(N_states,N_states,numel(s_y)+1);         % initialize for speed
    P_diag = zeros(N_states,numel(s_y)+1);                 % initialize for speed
    pred_error = zeros(1,numel(s_y));

    xi_hat(:,1) = xi_hat_init;                                % to ignore inital transient take the mean of the second half of the test data
%     xi_hat(9:end,1) = alpha_all;
    P_hat(:,:,1) = P_hat_init;


    anneal_on = 1;
    kappa_0 = 10000; % initial kappa
    t_end_anneal = N_samples/20; 
    
    y = s_y;

%         try
    for t=2:(length(s_y))
        xi_0p = squeeze(xi_hat(:,t-1));
        P_0p = squeeze(P_hat(:,:,t-1));
        % predict
        %N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q
        [xi_1m, P_1m] = prop_mean_and_cov(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q);
        if (t<=t_end_anneal) && anneal_on
            kappa = kappa_0^((t_end_anneal-t)/(t_end_anneal-1));
        else
            kappa = 1;
        end
        K = P_1m*H'/(H*P_1m*H' + kappa*R);
        % correct
        xi_hat(:,t) = xi_1m + K*(y(t-1) - H*xi_1m);
        P_hat(:,:,t) = (eye(N_states) - K*H)*P_1m;
        P_diag(:,t) = diag(squeeze(P_hat(:,:,t)));
        pred_error(t-1) = (y(t-1) - H*xi_1m); 
    end
%         catch me
%             display('something wrong in Ch'); 
%             continue;
%         end
    xi_hat = xi_hat(:,2:end);
    v_pyr_hat = xi_hat(1,:) + xi_hat(7,:) + xi_hat(9,:); % pyramidal membrane potential
    v_pyr_hat = double(v_pyr_hat/50);

%         v_es_hat = xi_hat(5,:); % excitatory stellate membrane potential
%         v_ii_hat = xi_hat(3,:); % inhibitory interneuron membrane potential
%     theta_hat = xi_hat(9:end,:);
%     %disp('here 4')
% %         Fs_final = 150;
% %         v_pyr_hat = single(resample(v_pyr_hat,Fs_final,Fs)); % resample from 400 to Fs_final Hz
% 
% %         v_es_hat = single(resample(v_es_hat,Fs_final,Fs));
% %         v_ii_hat = single(resample(v_ii_hat,Fs_final,Fs));
% 
%     input_hat = theta_hat(1,:);
%     aIP_hat = theta_hat(2,:);
%     aPI_hat = theta_hat(3,:);
%     aPE_hat = theta_hat(4,:);
%     aEP_hat = theta_hat(5,:);

%         save(['~/Local_data_folder/kalman_estimation/parameter_estimates_' num2str(ich) '.mat'],'v_pyr_hat','input_hat','aIP_hat','aPI_hat','aPE_hat','aEP_hat','-v7.3');
    
        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,C,N_states,N_syn,N_inputs,N_samples,xi,v0,varsigma,Q,R,H] = set_params(input,input_offset,fs,scaler,rr)
%function [A,B,C,N_states,N_syn,N_inputs,N_samples,v0,varsigma,Q,R,H,y,Fs] = set_params(input,input_offset)

% set_params.m

% set parameters - remove redundant states

scale = scaler;                         % this is to get states and derivatives on the same order of magnitude

% set units of membrane potentials (1 for mV, 0 for V)
%
mV = 1;
V2mVfactor = 1e3;

% time parameters
%
Fs = fs; % sampling rate per second
dt = 1/Fs;                                                                 % time step (s)
TimeOfSim = 50;                                                             % time of simulation (s)
N_samples = round(TimeOfSim/dt);                                            % number of time points to simulate

if ~isempty(input_offset)
    N_inputs = 2;
else
    N_inputs = 1; 
end

N_syn = 4;                            % number of synapses
N_states = 3*N_syn + N_inputs;        % u, z, alpha plus an input

%% define the disturbance covariance matrix
%
sigma_all = 5e-8;                               % something small for all states
sigma_input = 5e-4;                             % for input
sigma_params = 5e-5;%sigma_all;
sigma_offset = 5e-6;
sigma_R = 1e-3;

Q = eye(N_states)*(scale*sqrt(dt)*sigma_all)^2;                             % add a tiny bit of noise to all states (added for numerical stability)
Q(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*(scale*sqrt(dt)*sigma_params)^2;
% **** HARDCODED NUMBERS HERE
Q(2*N_syn+1,2*N_syn+1) = (scale*sqrt(dt)*sigma_input)^2;
if N_inputs > 1
    Q(2*N_syn+2,2*N_syn+2) = (scale*sqrt(dt)*sigma_offset)^2;
end

% measurement disturbance covariance
%
R = sigma_R^2;

%% General parameters from J&R
%
% sigmoid bits
%
f_max = 2.5;  % maximum firing rate (spikes/s)
r = rr; 
varsigma = 1.699/r;     % slope of the sigmoid                                                     % (spikes/(Vs))
varsigma_sq = varsigma^2;
v0 = 0.006;    % mean firing threshold                                                             % (V)

% synaptic gains
%
alpha_e = 3.25e-3;                                                          % gain of excitatory synapses (V)
alpha_i = -22e-3;                                                           % gain of inhibitory synapses (V)

% synaptic kernel time constants
%
ex_tau = 0.010;                     % excitatory synaptic time constant (s)
in_tau = 0.020;                     % inhibitory synaptic time constant (s)

% input to py population
%
input = input*scale;                   % ! scale here ! the input is a membrane potential input the pyramidal population this is similar to setting it at 270
input = input * alpha_e/ex_tau * ex_tau^2; % transformed input
%       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
%       input   synaptic gain    integral of kernel

% measurement DC offset
input_offset = input_offset * scale;        
input_offset = input_offset * alpha_e/ex_tau * ex_tau^2;


if mV == 1   % adjust units
    Q = V2mVfactor^2 * Q;
    R = V2mVfactor^2 * R;
    
    r = r/V2mVfactor;
    varsigma = 1.699/r;                                                     % (spikes/(Vs))
    varsigma_sq = varsigma^2;
    v0 = v0*V2mVfactor;
    alpha_e = alpha_e*V2mVfactor;                                           % gain of excitatory synapses (V)
    alpha_i = alpha_i*V2mVfactor;                                           % gain of inhibitory synapses (V)
    
    input= input*V2mVfactor;
    input_offset = input_offset*V2mVfactor;
end

% conectivity constants to relate to Jansen and Rit 1995 model
%
ConnectivityConst = 270;                            % Jansen and Rit connectivity parameters. Either 135, 270 or 675
C1 = ConnectivityConst;
C2 = 0.8*ConnectivityConst;
C3 = 0.25*ConnectivityConst;
C4 = 0.25*ConnectivityConst;

%% model structure
% ~~~~~~~~~~~~~~~~
%           X
%       __  |  __
%      /  \ | /  \
%     /  04 | 01  \
%     |     P     |
%  ^  |     | |   |  ^
%  |  E     | v   I  |  direction of info
%     03   /|\   02
%     |   / | \   |
%      \_/  |  \_/
%           v
% population types: E, P, I, X
% synapses: 01 (IP), 02 (PI), 03 (PE), 04 (EP)

%% this is the observation function.
%
H = zeros(1,N_states);        %Initialize to zeros and later add 1s to states that contribute to EEG

% initialize adjancy matrix
%
Gamma = zeros(2*N_syn+N_inputs,2*N_syn+N_inputs);   %  - plus 1 for input

% specify synapses
%
syn_index = 0;

% syn1, connection from I to P  % x = [ve ze vp1 zp1 vp2 zp2 vi zi];
%
syn_index = syn_index + 1;
tau(syn_index) = in_tau;
alpha(syn_index) = alpha_i*2*f_max*C4*dt / tau(syn_index);          % note the time constant and time step are in the gains
presyn_inputs = 2;                                                  % the presynaptic population is getting inputs from synapses 2
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;                       % set the entries of Gamma corresponding to indices of presynaptic inputs to 1
end
H(2*syn_index-1) = 1;

% syn2, connection from P to I
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
alpha(syn_index) = alpha_e*2*f_max*C3*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
presyn_inputs = [1 4 5];                                            % the presynaptic population is getting inputs from synapses 1, 4, 5
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 0;                                               % set to one if it contributes to the EEG (i.e. if the synapse is to Py cells)

% syn3, connection from P to E
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
alpha(syn_index) = alpha_e*2*f_max*C1*dt / tau(syn_index);
presyn_inputs = [1 4 5];                                         	% the presynaptic population is getting inputs from no other synapses (in the model)
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 0;

% syn4, connection from E to P
%
syn_index = syn_index + 1;
tau(syn_index) = ex_tau;
alpha(syn_index) = alpha_e*2*f_max*C2*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
presyn_inputs = 3;                                                  % the presynaptic population is getting inputs from synapse 3
if ~isempty(presyn_inputs)
    Gamma(2*syn_index,2*presyn_inputs-1) = 1;
end
H(2*syn_index-1) = 1;

% for input
%
syn_index = syn_index + 1;
H(2*syn_index-1) = 1;           % the input contributes to the observation function

if N_inputs > 1
    % offset term
    H(2*syn_index) = 1;            % offset contributes to the observation function
end

% rescale H
%
H = H/scale;                                                        % !scale! this is help deal with our numerical issues.

% set_ABC.m

% dean freestone

% this script take the parameters from set_params and creates the system
% matrices for the neural mass model (recurrent NN)

%% define A
% A encodes the dynamics induced by the membrane time constants.
% A is made up of the submatrices Psi in a block diagonal structure.
% There is a Psi for each connection in the model. This is where all the
% synaptic time constants enter the system. Further, the scale paramter
% enters here (and with C (multiplicative factor) and with the H (divisor).

Psi = zeros(2*N_syn,2*N_syn);               % initialise Psi, the component of A for fast states
for n=1:N_syn                               % build block diagonal structure
    index = 2*(n-1)+1;
    Psi(index:index+1,index:index+1) = [0 scale ; -1/(scale*tau(n)^2) -2/(tau(n))];
end

% A = [1 + dt*Psi, 0
%      0           1]
% where 1 is the identity matrix of the appropriate size
%
A = [(eye(2*N_syn) + dt*Psi) zeros(2*N_syn,N_syn+N_inputs) ; zeros(N_syn+N_inputs,2*N_syn) eye(N_syn+N_inputs)]; % top left part is discrete time versions.

%% define B (syanptic gain selection)
%
% Theta = [0 0 ... 0
%          1 0 ... 0
%          0 0 ... 0
%          0 1 ... 0
%          ...
%          0 0 ... 0
%          0 0 ... 1]
% B = [0 Theta ; 0 0]

Theta = zeros(2*N_syn,N_syn);                   % THETA IS USED TWICE
for n=1:N_syn
    index = 2*(n-1)+1;
    Theta(index:index+1,n) = [0 ; 1];
end
B = [zeros(2*N_syn,2*N_syn+N_inputs) Theta ; zeros(N_syn+N_inputs,3*N_syn+N_inputs)];

%% define C (adjacency matrix)
%
C = [Gamma/scale zeros(2*N_syn+N_inputs,N_syn) ; zeros(N_syn,3*N_syn+N_inputs)];

%%
%%
% WE ARE NOW BORN TO RUN
% ~~~~~~~~~~~~~~~~~~~~~~~
% set up the forward simulation
%

% define initial conditions
%
ve = 0; ze = 0; vp1 = 0; zp1 = 0; vp2 = 0; zp2 = 0; vp3 = 0; zp3 = 0; vi = 0; zi = 0;
x = [ve ze vp1 zp1 vp2 zp2 vi zi];

xi = zeros(N_states,N_samples);
xi(:,1) = [x input input_offset alpha]';                 % note: input and alpha are set in set_params

rng(1)
w = mvnrnd(zeros(1,N_states),Q,N_samples)';

phi_p = zeros(1,N_samples);

% run the model forward !!!This is the model
%
for n=1:N_samples-1
    
    phi = g(C*xi(:,n),v0,varsigma); % g is error function sigmoid
    phi_p(n) = phi(4);
    xi(:,n+1) = A*xi(:,n) + B*xi(:,n).*phi + w(:,n);
    
end

rng(1)
v = sqrt(R)*randn(1,N_samples); % variance
y = H*xi + v;                           % simulated EEG signal
disp(['std is ' num2str(std(y))])
end

function [xi_1m, P_1m] = prop_mean_and_cov(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q)

v_indexes = 1:2:2*N_syn+2;              % extra plus 2 for v_in the input!
z_indexes = 2:2:2*N_syn;
alpha_indexes = 2*N_syn+N_inputs+1:N_states;

w_fast = [0.1713244923791705 0.3607615730481384  0.4679139345726904 0.1713244923791705  0.3607615730481384 0.4679139345726904];
y_fast = [0.033765242898424 0.169395306766868 0.380690406958402 0.966234757101576 0.830604693233132 0.619309593041599];

CPC = C(z_indexes,:)*P_0p*C(z_indexes,:)';
dCPB = diag(C(z_indexes,:)*P_0p*B(z_indexes,:)');
Bxi = B(z_indexes,alpha_indexes)*xi_0p(alpha_indexes);
AP = A*P_0p;

%% analytic mean
% a priori state estimate
gamma = 1./sqrt(2*(diag(CPC) + varsigma^2));
beta = (C(z_indexes,v_indexes)*xi_0p(v_indexes) - v0).*gamma;  
try% <- v_0
    Xi = (erf(beta) + 1)/2;
catch
    disp(beta)
end
Upsilon = exp(-(beta.^2)).*gamma./sqrt(pi);
psi = Bxi.*Xi + dCPB.*Upsilon;                              % E[ Bxi_t o g(Cxi_t) ]

xi_1m = A*xi_0p;                                            % m is for minus (prior (or previous a posteriori)), p for plus (a priori)
xi_1m(2:2:2*N_syn) = xi_1m(2:2:2*N_syn) + psi;

%% exact covariance
%
% cov part 1 (Phi)
%
q2 = Upsilon.*(Bxi - dCPB.*beta.*gamma.^2 *2);
Phi = ones(N_states,1)*q2'.*(AP*C(z_indexes,:)') + ones(N_states,1)*Xi'.*(AP*B(z_indexes,:)');

% note: A term A*xi_t * E[ Bxi_t o g(Cxi_t) ]' = A*xi_0p*psi' is cancelled
% from here from the full covariance expansion.

% cov part 2 (Omega)
%
% NOTE: we can reduce the dimensionality (only compute upper triangle part) of the matrices we turn to vectors to increase speed in larger networks.
%
CPCgammagamma = asin(CPC.*(gamma*gamma')*2);
CPCgammagamma = CPCgammagamma(:);                                           % change to a vector
CPCgammagamma_y = CPCgammagamma*y_fast;                                     % a matrix of row vectors

betabeta = (beta*beta')*2;
betabeta_mat = betabeta(:)*ones(1,length(w_fast));                          % change to a vector and repmat

beta2mat = (beta.^2)*ones(1,N_syn);                                         % sq and repmat
beta2matT = beta2mat';                                                      % the transpose allow for permutations when we sum below
beta2_plus_beta2T = (beta2mat(:) + beta2matT(:))*ones(1,length(w_fast));    % change to a vectors, add, and repmat

% put it together
%
Psi = reshape(sum(CPCgammagamma*w_fast.*exp(-(beta2_plus_beta2T - betabeta_mat.*sin(CPCgammagamma_y))./cos(CPCgammagamma_y).^2),2),4,4)/(4*pi);
%                 ~~~~~~~~~~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%             ~~~
Omega = ((Xi*Xi') + Psi).*(Bxi*Bxi' + P_0p(alpha_indexes,alpha_indexes));
%        ~~~~~~~~~~~~~~    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         E[g^2(Cxi)]           (alpha^2 + sigma_alpha^2)

%% now allowing for correlations between params and states
% ~~~~~~~~
% v = sqrt(2*diag(d) + varsigma_sq);            % a strange + + was here: v= sqrt(2*diag(d) + + varsigma_sq); removed a plus
% phi = (erf(l.*varsigma./v) + 1)/2;
% z = e./f .* exp(-c.^2./(2*f)) .* ...
%     ((varsigma*e./(pi*v)) .* exp(1 - varsigma^3*c.^2./(2*sqrt(f).*v)) + phi.* (4*b.*f - 2*c.*e)./sqrt(2*pi*f));
% ~~~~~~~~
%%
% here we construct the cov mat.
%
P_1m = AP*A' + Q;
P_1m(z_indexes,z_indexes) = P_1m(z_indexes,z_indexes)  + Omega - psi*psi';
P_1m(:,z_indexes) = P_1m(:,z_indexes) + Phi;
P_1m(z_indexes,:) = P_1m(z_indexes,:) + Phi';


end

% error function sigmoid    
function out = g(v,v0,varsigma)

    out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
%       out = 1 ./ (1 + exp(varsigma*(-v+v0)));

end