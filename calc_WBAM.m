%% Calculate WBAM-related metrics and make all unitless

function [cycle,peak,int,inds] = calc_WBAM(WBAM,HS,TO,mass,height,plane)

    % Organize L and R HS's
    LHS_index = HS.L;
    RHS_index = HS.R;
    LTO_index = TO.L;
    RTO_index = TO.R;
    
    g = 9.81;
    norm = mass*(height.^(3/2))*sqrt(g);

    % LEFT
    L_WBAM_peak = nan.*ones(length(LHS_index),1);
    L_WBAM_int = nan.*ones(length(LHS_index),1);
    L_WBAM_cycle = nan.*ones(length(LHS_index)-1,100);
    L_WBAM_inds = nan.*ones(length(LHS_index),1);
    for ii = 1:length(LHS_index)
        % From HS to contralateral HS
        index1 = LHS_index(ii);
        index2 = RHS_index(find(RHS_index>LHS_index(ii),1,"first"));
        index3 = LTO_index(find(LTO_index>LHS_index(ii),1,"first"));

        % Peak
        if isempty(find(RHS_index>LHS_index(ii),1,"first"))
            L_WBAM_peak(ii,1) = nan;
            L_WBAM_inds(ii,1) = nan;
        else
            if contains(plane,'S')
                clearvars ind
                [L_WBAM_peak(ii,1), ind] = min((WBAM(index1:index2))./norm);
                L_WBAM_inds(ii,1) = ind/length(WBAM(index1:index2))*50; % percent of gait cycle
            elseif contains(plane,'F')
                clearvars ind
                [L_WBAM_peak(ii,1), ind] = max((WBAM(index1:index2))./norm);
                L_WBAM_inds(ii,1) = ind/length(WBAM(index1:index2))*50; % percent of gait cycle
            end
            if L_WBAM_peak(ii) == 0
                L_WBAM_peak(ii,1) = NaN;
                L_WBAM_inds(ii,1) = NaN;
            end
        end

        % Integrated
        if isempty(find(LTO_index>LHS_index(ii),1,"first"))
            L_WBAM_int(ii,1) = nan;
        else
            clearvars L_WBAM_step
            L_WBAM_step = (WBAM(index1:index3))./norm;
            L_WBAM_int(ii,1) = trapz(L_WBAM_step)/100; % divide by 100 Hz or 1/s and we need to integrate over time, not frames
        end
    end

    for ii = 1:length(LHS_index)-1
        index1 = LHS_index(ii);
        index2 = LHS_index(find(LHS_index>LHS_index(ii),1,"first"));
        L_WBAM_stride{ii,1} = (WBAM(index1:index2))./norm;
    end
    L_WBAM_cycle(:,1:100) = resample_gait_cycle(L_WBAM_stride);

    % RIGHT
    R_WBAM_peak = nan.*ones(length(RHS_index),1);
    R_WBAM_int = nan.*ones(length(RHS_index),1);
    R_WBAM_cycle = nan.*ones(length(RHS_index)-1,100);
    R_WBAM_inds = nan.*ones(length(RHS_index),1);
    for ii = 1:length(RHS_index)
        % From HS to contralateral HS
        index1 = RHS_index(ii);
        index2 = LHS_index(find(LHS_index>RHS_index(ii),1,"first"));
        index3 = RTO_index(find(RTO_index>RHS_index(ii),1,"first"));

        % Peak
        if isempty(find(LHS_index>RHS_index(ii),1,"first"))
            R_WBAM_peak(ii,1) = nan;
            R_WBAM_inds(ii,1) = nan;
        else
            clearvars ind
            [R_WBAM_peak(ii,1), ind] = min((WBAM(index1:index2))./norm);
            R_WBAM_inds(ii,1) = ind/length(WBAM(index1:index2))*50; % percent of gait cycle
            if R_WBAM_peak(ii) == 0
                R_WBAM_peak(ii,1) = NaN;
                R_WBAM_inds(ii,1) = NaN;
            end
        end

        % Integrated
        if isempty(find(RTO_index>RHS_index(ii),1,"first"))
            R_WBAM_int(ii,1) = nan;
        else
            clearvars R_WBAM_step
            R_WBAM_step = (WBAM(index1:index3))./norm;
            R_WBAM_int(ii,1) = trapz(R_WBAM_step)/100; % divide by 100 Hz or 1/s and we need to integrate over time, not frames
        end
    end

    for ii = 1:length(RHS_index)-1
        index1 = RHS_index(ii);
        index2 = RHS_index(find(RHS_index>RHS_index(ii),1,"first"));
        R_WBAM_stride{ii,1} = (WBAM(index1:index2))./norm;
    end
    R_WBAM_cycle(:,1:100) = resample_gait_cycle(R_WBAM_stride);

    % Sort temporally
    [~,sort_ind] = sort([LHS_index; RHS_index]);

    temp_WBAM_peak = [L_WBAM_peak; R_WBAM_peak];
    peak.sorted(:,1) = temp_WBAM_peak(sort_ind);
    peak.L = L_WBAM_peak;
    peak.R = R_WBAM_peak;

    temp_WBAM_int = [L_WBAM_int; R_WBAM_int];
    int.sorted(:,1) = temp_WBAM_int(sort_ind);
    int.L = L_WBAM_int;
    int.R = R_WBAM_int;

    cycle.L = L_WBAM_cycle;
    cycle.R = R_WBAM_cycle;

    inds.L = L_WBAM_inds;
    inds.R = R_WBAM_inds;

end