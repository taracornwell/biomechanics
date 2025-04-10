%% Function for step width, length, time

function [SW, SL, ST] = spatiotemp(HS, MarkerPos)

    % Deconstruct the structures
    LHS_index = HS.L;
    RHS_index = HS.R;
    L_pos = MarkerPos.L;
    R_pos = MarkerPos.R;
    
    % Right or left step first?
    if LHS_index(1) > RHS_index(1) % right step first
        first = 'R';
    else
        first = 'L';
    end

    % LEFT
    for ii = 1:length(LHS_index)
        if ii == 1 && first == 'L'
            % skip because no previous step
            SW.L(ii,1) = nan;
            SL.L(ii,1) = nan;
            ST.L(ii,1) = nan;
        else
            ind2 = LHS_index(ii); % current LHS
            ind1 = RHS_index(find(RHS_index < ind2,1,'last')); % last RHS that occurred before current LHS
            
            % Step width = ML distance between markers at HSs
            SW.L(ii,1) = L_pos(ind2,2) - R_pos(ind1,2); % y is ML and + to the left

            % Step length = AP distance between markers at current HS
            SL.L(ii,1) = L_pos(ind2,1) - R_pos(ind2,1); % x is AP and + forward

            % Step time = time (s) between HSs
            ST.L(ii,1) = (double((ind2 - ind1)))./100;
        end
    end

    % RIGHT
    for ii = 1:length(RHS_index)
        if ii == 1 && first == 'R'
            % skip because no previous step
            SW.R(ii,1) = nan;
            SL.R(ii,1) = nan;
            ST.R(ii,1) = nan;
        else
            ind2 = RHS_index(ii); % current RHS
            ind1 = LHS_index(find(LHS_index < ind2,1,'last')); % last LHS that occurred before current RHS
            
            % Step width = ML distance between markers at HSs
            SW.R(ii,1) = L_pos(ind1,2) - R_pos(ind2,2); % y is ML and + to the left

            % Step length = AP distance between markers at current HS
            SL.R(ii,1) = R_pos(ind2,1) - L_pos(ind2,1); % x is AP and + forward

            % Step time = time (s) between HSs
            ST.R(ii,1) = (double((ind2 - ind1)))./100;
        end
    end

    % Sort temporally
    [~,sort_ind] = sort([LHS_index; RHS_index]);

    temp_SW = [SW.L; SW.R];
    SW.sorted(:,1) = temp_SW(sort_ind);

    temp_SL = [SL.L; SL.R];
    SL.sorted(:,1) = temp_SL(sort_ind);

    temp_ST = [ST.L; ST.R];
    ST.sorted(:,1) = temp_ST(sort_ind);

end