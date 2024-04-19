function eq = find_equilibrium_network(x0,options,n_channels,aIP_arr,aPI_arr,aPE_arr,aEP_arr,w_matrix,tau_e,tau_i,H,tau_d,v0,varsigma,scale)
    fun = @system;
    try
        eq = fsolve(fun,x0,options);
    catch
        x0 = zeros(1, n_channels*10);
        eq = fsolve(fun,x0,options);
    end
    
    function F = system(x)
        for i = 1:n_channels
            j = 10*(i-1);
            
            % define model states and parameters
            vip = x(j+1);
            zip = x(j+2);
            vpi = x(j+3);
            zpi = x(j+4);
            vpe = x(j+5);
            zpe = x(j+6);
            vep = x(j+7);
            zep = x(j+8);
            mu = x(j+9);
            u = x(j+10);
            aIP = aIP_arr(i);
            aPI = aPI_arr(i);
            aPE = aPE_arr(i);
            aEP = aEP_arr(i);
            w = w_matrix(:,i);
            
            % define model ode
            F(j+1) = scale*zip;
            F(j+2) = aIP*g(vpi/scale) - 2/tau_i*zip - vip/(scale*tau_i^2);
            
            F(j+3) = scale*zpi;
            F(j+4) = aPI*g((vip+vep+mu)/scale) - 2/tau_e*zpi - vpi/(scale*tau_e^2);
            
            F(j+5) = scale*zpe;
            F(j+6) = aPE*g((vip+vep+mu)/scale) - 2/tau_e*zpe - vpe/(scale*tau_e^2);
            
            F(j+7) = scale*zep;
            F(j+8) = aEP*g(vpe/scale) - 2/tau_e*zep - vep/(scale*tau_e^2);
            
            F(j+9) = scale*u;
            F(j+10) = g(x*H/scale)*w - 2/tau_d*u - mu/(scale*tau_d^2);
        end
    end

    function out = g(v)
        out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
    end
end