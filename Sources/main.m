function main(tasks)
    if ismember("parameters estimation", tasks) % run whole-cortex model parameters estimation
        disp('Running the 1st showcase: whole-cortex model parameters estimation...')
        tic
        results = evalc('calculate_dynamic_stability');
        toc
    end  
    
    if ismember("cortical stability", tasks) % run dynamic cortical stability
        disp('Running the 2nd showcase: dynamic cortical stability...')
        tic
        % results = evalc('calculate_dynamic_stability');
        results = evalc('show_dynamic_stability');
        toc
    end
end
