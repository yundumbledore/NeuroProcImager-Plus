function J = construct_Jacobian(n_channels, eq, w_matrix, aIP_array, aPI_array, aPE_array, aEP_array, tau_e, tau_i, tau_d, varsigma, v0, scale)
    J = zeros(n_channels*10, n_channels*10);
    
    for i = 1:n_channels
        aIP = aIP_array(i);
        aPI = aPI_array(i);
        aPE = aPE_array(i);
        aEP = aEP_array(i);
        w = w_matrix(:,i);
        
        for j = 1:n_channels
            vip = eq((j-1)*10+1);
            zip = eq((j-1)*10+2);
            vpi = eq((j-1)*10+3);
            zpi = eq((j-1)*10+4);
            vpe = eq((j-1)*10+5);
            zpe = eq((j-1)*10+6);
            vep = eq((j-1)*10+7);
            zep = eq((j-1)*10+8);
            mu = eq((j-1)*10+9);
            u = eq((j-1)*10+10);
            
            if i == j
                sigma1 = -1/(scale*tau_e^2);
                sigma5 = exp(-(0.5*(v0-(mu+vep+vip)/(scale))^2)/(varsigma^2));
                sigma2 = 0.3989*w(j)*sigma5/(scale*varsigma);
                sigma3 = 0.3989*aPI*sigma5/(scale*varsigma);
                sigma4 = 0.3989*aPE*sigma5/(scale*varsigma);
                
                
                J((i-1)*10+1:i*10, (j-1)*10+1:j*10) = [0 scale 0 0 0 0 0 0 0 0;
                    -1/(scale*tau_i^2) -2/tau_i (0.3989*aIP*exp((-0.5*(v0-(vpi)/(scale))^2)/(varsigma^2)))/(scale*varsigma) 0 0 0 0 0 0 0;
                    0 0 0 scale 0 0 0 0 0 0;
                    sigma3 0 sigma1 -2/tau_e 0 0 sigma3 0 sigma3 0;
                    0 0 0 0 0 scale 0 0 0 0;
                    sigma4 0 0 0 sigma1 -2/tau_e sigma4 0 sigma4 0;
                    0 0 0 0 0 0 0 scale 0 0;
                    0 0 0 0 (0.3989*aEP*exp((-0.5*(v0-(vpe)/(scale))^2)/(varsigma^2)))/(scale*varsigma) 0 sigma1 -2/tau_e 0 0;
                    0 0 0 0 0 0 0 0 0 scale;
                    sigma2 0 0 0 0 0 sigma2 0 sigma2-1/(scale*tau_d^2) -2/tau_d];
            else
                sigma5 = exp(-(0.5*(v0-(mu+vep+vip)/(scale))^2)/(varsigma^2));
                sigma2 = 0.3989*w(j)*sigma5/(scale*varsigma);
                J((i-1)*10+1:i*10, (j-1)*10+1:j*10) = [zeros(9,10);
                    sigma2 0 0 0 0 0 sigma2 0 sigma2 0];
            end
        end
    end
end