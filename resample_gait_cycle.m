% Input signal should be cells
function resampled = resample_gait_cycle(signal)

    resampled = nan.*ones(length(signal),100);
    
    for ii = 1:length(signal)
        clearvars flipped resamp_long
        flipped = [flip(signal{ii});signal{ii};flip(signal{ii})];
        resamp_long = resample(flipped,300,length(flipped));
        resampled(ii,:) = resamp_long(101:200)';
    end
end