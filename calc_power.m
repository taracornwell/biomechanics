%% Calculate mechanical power done by leg per step
% Tara Cornwell - Summer 2023
% function [L_power,R_power,L_power_bod,L_power_TM,R_power_bod,R_power_TM] 
%           = calc_power(LHS, LTO, RHS, RTO, Lforce, Rforce, VCOM, TMspeeds, pert_events)
% Inputs:   
%           LHS, LTO, RHS, RTO with indices of heel-strike and toe-off events
%           Lforce and Rforce with force data from current trial (N, sampled at 1000Hz)
%           VCOM = COM velocity from current trial in 3 columns (m/s)
%           TMspeeds has 2 structures SSWS, pert w/ treadmill speeds (m/s)
%           pert_events has 2 structures L, R w/ LHS/RHS indices denoting
%               which steps are the perturbation steps
% Outputs:
%           L_power & R_power have a cell per step with instantaneous
%               values of power per frame (W)
%           L/R_power_bod and _TM have a cell per step with instantaneous
%           values of the leg's power on the body and treadmill

function [L_power,R_power,L_power_bod,L_power_TM,R_power_bod,R_power_TM] = calc_power(LHS,...
    LTO, RHS, RTO, Lforce, Rforce, VCOM, TMspeeds, pert_events)

    % Downsample first bc sampled at 1kHz, not 100Hz

    div = length(Lforce)/length(VCOM);
    Lforce_ds = downsample(Lforce,div);
    Rforce_ds = downsample(Rforce,div);

    % Loop through left steps
    for ii = 1:length(LHS)
        clearvars start finish F V_bod V_TM

        start = LHS(ii);                                % index for current heel-strike
        finish = LTO(find(LTO>LHS(ii),1,"first"));  % index for corresponding toe-off

        % Check if current step is perturbation step and set treadmill
        % speed accordingly
        if ismember(ii,pert_events.L)
            V_TM(1,1:3) = [-TMspeeds.pert, 0, 0];
        else
            V_TM(1,1:3) = [-TMspeeds.SSWS, 0, 0];
        end

        % Pull force and velocity values for just the current step
        F = Lforce_ds(start:finish, 1:3);
        V_bod = VCOM(start:finish, 1:3);
    
        % Loop through every frame
        for jj = 1:size(F,1)
    
            % Calculate instantaneous leg power at every frame
            % Power of leg on the body
            L_power_bod{ii,1}(jj,1) = dot(F(jj,1:3), V_bod(jj,1:3));
            % Power of leg on treadmill
            L_power_TM{ii,1}(jj,1) = dot(-F(jj,1:3), V_TM(1,1:3));
    
            % Sum of instantaneous powers (power of leg on body + leg on tread) = TOTAL POWER
            L_power{ii,1}(jj,1) = L_power_bod{ii}(jj) + L_power_TM{ii}(jj);

        end

    end

    % Loop through left steps
    for ii = 1:length(RHS)
        clearvars start finish F V_bod V_TM

        start = RHS(ii);                                % index for current heel-strike
        finish = RTO(find(RTO>RHS(ii),1,"first"));  % index for corresponding toe-off

        % Check if current step is perturbation step and set treadmill
        % speed accordingly
        if ismember(ii,pert_events.R)
            V_TM(1,1:3) = [-TMspeeds.pert, 0, 0];
        else
            V_TM(1,1:3) = [-TMspeeds.SSWS, 0, 0];
        end

        % Pull force and velocity values for just the current step
        F = Rforce_ds(start:finish, 1:3);
        V_bod = VCOM(start:finish, 1:3);
    
        % Loop through every frame
        for jj = 1:size(F,1)
    
            % Calculate instantaneous leg power at every frame
            % Power of leg on the body
            R_power_bod{ii,1}(jj,1) = dot(F(jj,1:3), V_bod(jj,1:3));
            % Power of leg on treadmill
            R_power_TM{ii,1}(jj,1) = dot(-F(jj,1:3), V_TM(1,1:3));
    
            % Sum of instantaneous powers (power of leg on body + leg on treadmill) = TOTAL POWER
            R_power{ii,1}(jj,1) = R_power_bod{ii}(jj) + R_power_TM{ii}(jj);

        end
        
    end

end