function cycles = separate_by_cycle(signal,HS,left,right)

    % Organize L and R HS's
    LHS_index = HS.L;
    RHS_index = HS.R;

    % LEFT
    for ii = 1:length(LHS_index)-1
        % From HS to TO (ipsilateral)
        index1 = LHS_index(ii);
        index2 = LHS_index(find(LHS_index>LHS_index(ii),1,"first"));

        if isempty(find(LHS_index>LHS_index(ii),1,"first"))
            L_cycle{ii,1} = nan; % to make same length as other vectors
        else
            L_cycle{ii,1} = signal(index1:index2);
        end
    end

    % RIGHT
    for ii = 1:length(RHS_index)-1

        % From HS to TO (ipsilateral)
        index1 = RHS_index(ii);
        index2 = RHS_index(find(RHS_index>RHS_index(ii),1,"first"));

        if isempty(find(RHS_index>RHS_index(ii),1,"first"))
            R_cycle{ii,1} = nan; % to make same length as other vectors
        else
            R_cycle{ii,1} = signal(index1:index2);
        end
    end

    if left == 1
        cycles.L(:,1:100) = resample_gait_cycle(L_cycle);
    end
    if right == 1
        cycles.R(:,1:100) = resample_gait_cycle(R_cycle);
    end

end