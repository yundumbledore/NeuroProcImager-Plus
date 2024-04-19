% NMM_DEFINE Neural Mass Model based on Jensen and Rit. Defines the state based representation of the model.
%
% Inputs:
%   x0 - initial states
%   P0 - initial covariance matrix
%   params - parameters defined by Jensen and Rit
%
% Outputs:
%   nmm - (struct) the neural mass model
%
function nmm = nmm_define(x0,P0,params)
% Indexes
v_idx = [1 3 5 7];
z_idx = [2 4 6 8];
u_idx = 9;
alpha_idx = [10 11 12 13];

% the parameters
dt          = params.dt;
e_0         = params.e0;
r           = params.r;	% varsigma
v0          = params.v0; % Threshold
decay_e     = params.decay_e; % inverse time constants (excitatory)
decay_i     = params.decay_i; % (inhibitory)
alpha_ei    = params.alpha_ei; % synaptic gains (excitatory)
alpha_ie    = params.alpha_ie; % (inhibitory)
u           = params.u;	% mean input firing rate.
scale       = params.scale; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2

c_constant = 270; %100;
c1 = 1*c_constant;	% number of synapses
c2 = 0.8*c_constant;
c3 = 0.25*c_constant;
c4 = 0.25*c_constant;

% Number of augmented states
xlen = length(x0);

% Linear component of model
A =     [1,                  dt*scale,      0,              0,          0,  0,  0, 0, 0, 0, 0, 0, 0; ...
  -decay_i^2*dt/scale,	1-2*decay_i*dt,     0,              0,          0,  0,  0, 0, 0, 0, 0, 0, 0; ...
         0,                     0,          1,            dt*scale,     0,  0,  0, 0, 0, 0, 0, 0, 0; ...
         0,                     0  -decay_e^2*dt/scale,	1-2*decay_e*dt,	0,  0,  0, 0, 0, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          1,  dt*scale,  0, 0, 0, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          -decay_e^2*dt/scale,  1-2*decay_e*dt,  0, 0, 0, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  1, dt*scale, 0, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  -decay_e^2*dt/scale, 1-2*decay_e*dt, 0, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  0, 0, 1, 0, 0, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  0, 0, 0, 1, 0, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  0, 0, 0, 0, 1, 0, 0; ...
         0,                     0,          0,              0,          0,  0,  0, 0, 0, 0, 0, 1, 0; ...
         0,                     0,          0,              0,          0,  0,  0, 0, 0, 0, 0, 0, 1]; ...
     
% B Matrix (Augmented parameters)                                                                
%      
% B =     [0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	1,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	1; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,  0,  0];
%
B = zeros(xlen);
B(z_idx, alpha_idx) = diag(ones(size(z_idx)));

% C Matrix (Augmented)
% C =     [0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	1,	0,	1,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          1,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,  0,  0];
%
C = zeros(xlen);
C(2, 3) = 1; 
C(4, [1 7 9]) = 1; 
C(6, [1 7 9]) = 1;
C(8, 5) = 1;
C = C./scale;

alpha_ip = alpha_ie * c4 * 2 * e_0 * dt * decay_i; %
alpha_pi = alpha_ei * c3 * 2 * e_0 * dt * decay_e; % 
alpha_pe = alpha_ei * c1 * 2 * e_0 * dt * decay_e; % 
alpha_ep = alpha_ei * c2 * 2 * e_0 * dt * decay_e; % 

% SCALE 1 - this is to avoid large differences between states upsetting the filter 
% (magnitude of the membrane potentials and their derivatives)
input = scale*u;
% SCALE 2 - this converts a constant input to its effect on the pyramidal
% membrane potential by taking the steady state limit of the synaptic kernel
% (assumption that the input varies much slower than the state variables).
% input = input * alpha_ei*decay_e /decay_e^2;
%       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
%       input   synaptic gain    integral of kernel

x0(9) = input;
x0(10) = alpha_ip;
x0(11) = alpha_pi;
x0(12) = alpha_pe;
x0(13) = alpha_ep;

nmm = struct;
nmm.A = A;
nmm.B = B;
nmm.C = C;
nmm.x0 = x0;
nmm.P0 = P0;
nmm.params = params;
nmm.options = struct;

end % end function - nmm_define
