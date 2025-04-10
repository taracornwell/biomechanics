%% Get indices of HS and TO events

% Inputs: 4 vectors of time-series data of the HS, TO events for left and
% right sides
% Outputs: 2 structures of HS and TO index data organized by left (.L) and 
% right (.R) sides

function [HS, TO] = get_event_inds(LHS, LTO, RHS, RTO)

        % Get indices
        LHS_index(:,1) = int32(LHS*100+1); % index = time*freq +1
        RHS_index(:,1) = int32(RHS*100+1);
        LTO_index(:,1) = int32(LTO*100+1);
        RTO_index(:,1) = int32(RTO*100+1);

        % Must start with HS and end with TO
        if LTO_index(1,1) < LHS_index(1,1) % if first TO event happens before first HS event
            LTO_index = LTO_index(2:end); % remove it
        end
        if RTO_index(1,1) < RHS_index(1,1)
            RTO_index = RTO_index(2:end);
        end
        if LHS_index(end,1) > LTO_index(end,1) % if last HS event happens after last TO event
            LHS_index = LHS_index(1:end-1,1); % remove it
        end
        if RHS_index(end,1) > RTO_index(end,1)
            RHS_index = RHS_index(1:end-1,1);
        end
        
        % Organize output variables
        HS.L = LHS_index;
        HS.R = RHS_index;
        TO.L = LTO_index;
        TO.R = RTO_index;

end