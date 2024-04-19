for sub_ind = 7:28
    for c = [3 7]
        try
            load(['~/cn25_scratch/_variable_estimates/variable_estimates_errorReport_' num2str(sub_ind) '_2_' num2str(c) '.mat'], 'chan_error')
            s = sum(chan_error);
            if s > 0
                disp([num2str(sub_ind) '_2_' num2str(c) ': ' num2str(s) '...'])
            end
        catch
            continue
        end
    end
end