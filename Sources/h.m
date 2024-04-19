function out = h(t)
    out = heaviside(t)*(t/(1/33))*exp((-t)/(1/33));
end