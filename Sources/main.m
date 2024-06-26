function main(tasks, file_name)
    if ismember("cortical stability", tasks) % run dynamic cortical stability
        
        disp('Running the showcase: dynamic cortical stability...')
        tic
        results = evalc('calculate_dynamic_stability(file_name)');
        results = evalc('show_dynamic_stability');
        toc
    end
end
