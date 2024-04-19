% Set parameters for the JR model
%
function params = set_parameters(mode, varargin)

mu = 300; % Default external input value

if nargin < 1
    mode = 'alpha';
elseif nargin > 1
    % If the function is called with 2 input arguments, the input (u) is 
    % defined by the second argument
    mu = varargin{1};
end

switch mode
    case 'alpha' % 8 Hz - 12 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope % params.r = 0.56;  % logistic sigmoid slope
        params.v0 = 6; % Firing Threshold

        % inverse time constants
        params.decay_e = 100;% (1/ tau_e)
        params.decay_i = 50; % (1/tau_i)
        
        params.alpha_ei = 3.25;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = -22;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000

        params.u = mu;%11;%220;%15;%11;        % mean input mem potential

        params.dt = 0.0025;     % sampling time step   
        
        params.scale = 50; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'beta' % 12 Hz - 30 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 150;% (1/ tau_e)
        params.decay_i = 150; % (1/tau_i)        
        params.alpha_ei = 22;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 12;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step           
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'gamma' % > 30 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 150*2;% (1/ tau_e)
        params.decay_i = 150*4; % (1/tau_i)
        params.alpha_ei = 22;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 12;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'delta' % 1 Hz - 4 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 100/4;% (1/ tau_e)
        params.decay_i = 50/3; % (1/tau_i)
        params.alpha_ei = 3.25;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 22;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
    otherwise
        error('%s rythm not implemented, sorry!', mode);
end
