%% Calculate mechanical work done by leg per step
% Tara Cornwell - Summer 2023
% function [L_pos_work,L_neg_work,R_pos_work,R_neg_work] = calc_work(L_power, R_power)
% Inputs:   
%           L_power & R_power have a cell per step with instantaneous
%               values of power per frame (W)
% Outputs:
%           L_pos_work & L_neg_work are vectors w/ + and - work values 
%               for the left leg per step - one value per left step (J)
%           R_pos_work & R_neg_work are vectors w/ + and - work values 
%               for the right leg per step - one value per right step (J)

function [L_pos_work,L_neg_work,R_pos_work,R_neg_work] = calc_work(L_power, R_power, Hz)

    if nargin < 3
        Hz = 100;
    end

    % Loop through left steps
    for ii = 1:length(L_power)

        % Integrate power to get + and - work performed per step
        posVector = [];
        negVector = [];

        %Figure out # of zero crossings to define pos and neg areas
        clearvars signs crossings num_crossings
        signs = sign(L_power{ii});
        crossings = diff(sign(signs))/2; % need to divide by 2 because -1 to +1 = diff of 2 but that's just 1 sign change
        num_crossings = numel(find(crossings));

        zero_cross_ind = 1; % this is 1 because area 1 is bounded by the first value and the first crossing
        for area = 1:num_crossings + 1 % loop through number of areas which is the # of crossings + 1
            clearvars ind1 next_cross
            ind1 = zero_cross_ind(area);

            if L_power{ii}(ind1) > 0  % signal starts positive
                next_cross = find(L_power{ii}(ind1:end) < 0, 1, 'first')-1+ind1; % find index where the signal is last positive. This bounds the current area
                if isempty(next_cross)
                    next_cross = length(L_power{ii}); % end of last area will probably not cross zero
                end
                % Calculate the current area = positive work
                posVector = [posVector; trapz(L_power{ii}(ind1:next_cross))];
            else
                % Current area starts negative
                next_cross = find(L_power{ii}(ind1:end) > 0, 1, 'first')-1+ind1; % find index where the signal is last negative. This bounds the current area
                if isempty(next_cross)
                    next_cross = length(L_power{ii}); % end of last area will probably not cross zero
                end
                % Calculate the current area = positive work
                negVector = [negVector; trapz(L_power{ii}(ind1:next_cross))];
            end
            zero_cross_ind = [zero_cross_ind; next_cross];
        end
         
        % Calculate positive and negative work for the current step (ii)
        L_pos_work(ii,1) = (sum(posVector,'omitnan'))/Hz; % get into J by dividing by sample rate 100 Hz
        L_neg_work(ii,1) = (-sum(negVector,'omitnan'))/Hz; % make positive so magnitude of neg work

    end

    % Loop through left steps
    for ii = 1:length(R_power)
        
        % Integrate power to get + and - work performed per step
        posVector = [];
        negVector = [];

        %Figure out # of zero crossings to define pos and neg areas
        clearvars signs crossings num_crossings
        signs = sign(R_power{ii});
        crossings = diff(sign(signs))/2; % need to divide by 2 because -1 to +1 = diff of 2 but that's just 1 sign change
        num_crossings = numel(find(crossings));

        zero_cross_ind = 1; % this is 1 because area 1 is bounded by the first value and the first crossing
        for area = 1:num_crossings + 1 % loop through number of areas which is the # of crossings + 1
            clearvars ind1 next_cross
            ind1 = zero_cross_ind(area);

            if R_power{ii}(ind1) > 0  % signal starts positive
                next_cross = find(R_power{ii}(ind1:end) < 0, 1, 'first')-1+ind1; % find index where the signal is last positive. This bounds the current area
                if isempty(next_cross)
                    next_cross = length(R_power{ii}); % end of last area will probably not cross zero
                end
                % Calculate the current area = positive work
                posVector = [posVector; trapz(R_power{ii}(ind1:next_cross))];
            else
                % Current area starts negative
                next_cross = find(R_power{ii}(ind1:end) > 0, 1, 'first')-1+ind1; % find index where the signal is last negative. This bounds the current area
                if isempty(next_cross)
                    next_cross = length(R_power{ii}); % end of last area will probably not cross zero
                end
                % Calculate the current area = positive work
                negVector = [negVector; trapz(R_power{ii}(ind1:next_cross))];
            end
            zero_cross_ind = [zero_cross_ind; next_cross];
        end
         
        % Calculate positive and negative work for the current step (ii)
        R_pos_work(ii,1) = (sum(posVector,'omitnan'))/Hz; % get into J by dividing by sample rate 100 Hz
        R_neg_work(ii,1) = (-sum(negVector,'omitnan'))/Hz; % make positive so magnitude of neg work
    
    end

end