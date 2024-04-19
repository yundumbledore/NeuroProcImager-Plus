 function [xi_hat, er] = AKF_quick(y, xi0, P0, A, B, C, Q, R, H, varsigma, v0, N_syn, N_states, N_inputs)
    xi_hat = zeros(N_states,numel(y)+1);
    P_hat = zeros(N_states,N_states,numel(y)+1); 

    xi_hat(:,1) = xi0;
    P_hat(:,:,1) = P0;

    anneal_on = 1;
    kappa_0 = 10000; % initial kappa
    t_end_anneal = 20000/20; 
    
    try
        for t=2:(length(y))
            xi_0p = squeeze(xi_hat(:,t-1));
            P_0p = squeeze(P_hat(:,:,t-1));
            % predict
            [xi_1m, P_1m] = prop_mean_and_cov(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q);
            if (t<=t_end_anneal) && anneal_on
                kappa = kappa_0^((t_end_anneal-t)/(t_end_anneal-1));
            else
                kappa = 1;
            end
            K = P_1m*H'/(H*P_1m*H' + kappa*R);
            % correct
            xi_hat(:,t) = xi_1m + K*(y(t-1) - H*xi_1m);
            % prevent IPSP amplitude from being positive
            if xi_hat(10,t) > 0
                xi_hat(10,t) = 0;
            end

            % prevent EPSP amplitude from being negative
%             if xi_hat(9,t) < 0
%                 xi_hat(9,t) = 0;
%             end

            if xi_hat(11,t) < 0
                xi_hat(11,t) = 0;
            end

            if xi_hat(12,t) < 0
                xi_hat(12,t) = 0;
            end

            if xi_hat(13,t) < 0
                xi_hat(13,t) = 0;
            end
            
            P_hat_corrected = (eye(N_states) - K*H)*P_1m;
%            % Force symmetry on P_hat
%             P_hat_corrected = (P_hat_corrected + P_hat_corrected')/2;
%             % Check eigenvalues
%             [~,flag] = chol(P_hat_corrected);
%             if flag
%                 % If any is negative, find the nearest Semipositve Definite
%                 % matrix
%                 try
%                     [P_hat_corrected, k]= nearestSPD(P_hat_corrected); % Nearest SPD, Higham (1988) - Parvin's method
%                 catch ME
%                     if strcmp('MATLAB:svd:matrixWithNaNInf', ME.identifier)
%                         disp(['Error during iteration: ?']);% num2str(n)]);
%                     end
%                     rethrow(ME);
%                 end
%                 if k == -1
%                     % Infinite loop in the nearestSPD script. No SPD matrix
%                     % found
%                     %disp(['Couldn''t find nearest SPD']);% at t = ?']);% num2str(n)]);
%                     error('Couldn''t find the nearest SPD');
%                 end
%             end
            P_hat(:,:,t) = P_hat_corrected;
        end
        
        er = 0;
    catch
        % disp(ME)
        er = 1;
    end
end

%% SUBFUNCTIONS
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
Xi = (erf(beta) + 1)/2;
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