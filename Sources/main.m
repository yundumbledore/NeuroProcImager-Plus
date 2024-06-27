function main(tasks)
    if ismember("cortical stability", tasks) % run dynamic cortical stability
        disp('Running the 2nd showcase: dynamic cortical stability...')
        tic
        results = evalc('calculate_dynamic_stability');
        results = evalc('show_dynamic_stability');
        toc
    end
end
