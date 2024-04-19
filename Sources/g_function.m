%% sigmoid function
function [out] = g_function(v,v0,varsigma)
    out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
end