function B = analytic_multiregress(X, Y)
    B = (X'*X)\(X')*Y;
    B(isnan(B)) = 0;
end