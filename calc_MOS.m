%% Margins of stability Function
% Written by Tara Cornwell Spring 2022
% MOS = BOS - XCOM
% Both anterior and lateral margins
% INPUTS:   HS and TO event indices within L and R structures,
%           Anterior and Lateral foot positions(:,3)
%           COM(:,3) and VCOM(:,3)
%           L0 = leg length (m)
%           TMspeeds structures: SSWS and pert (m/s),
%           pert_events structures: 
%               If perturbations exist: L and R with index numbers
%               indicating with L and R HS events are perturbation steps
%               If no perturbations: []
% OUTPUTS:  All output structures are separated by L and R
%           MOS structures: ant, lat with further structures: min, HS
%           BOS structures: ant, lat
%           XCOM structures: ant, lat
%           COM wrt BOS with structures: ant, lat

function [MOS, BOS, XCOM, COM, VCOM] = calc_MOS(HS,TO,AntFoot,LatFoot,COMpos,COMvel,L0,TMspeeds,pert_events)
    
    % Deconstruct the structures
    LHS_index = HS.L;
    LTO_index = TO.L;
    RHS_index = HS.R;
    RTO_index = TO.R;

    % Set variables depending on if perturbations occur
        if isempty(pert_events) == 0 % perturbations occur
            perts = 1; % if perturbations occur = 1; otherwise = 0
        else % no perturbations labeled
            perts = 0;
        end

    % Define constants
    g = 9.81; % gravitational acceleration (m/s)

    % Angular frequency
    w = sqrt(g/L0);
    
    % Define BOS and XCOM = COM + VCOM/w at each HS
    % LEFT
    for ii = 1:length(LHS_index)
        % From HS to TO (ipsilateral)
        ind1 = LHS_index(ii);
        ind2 = LTO_index(find(LTO_index>ind1,1,'first'));

        % Assign BOS based on anterior/lateral foot positions
        BOS.ant.L{ii,1} = AntFoot.L(ind1:ind2,1); % X is AP
        BOS.lat.L{ii,1} = LatFoot.L(ind1:ind2,2); % Y is ML
        
        % Calculate XCOM = COM + VCOM/w
        XCOM.lat.L{ii,1} = COMpos(ind1:ind2,2) + COMvel(ind1:ind2,2)./w;
        COM.lat.L{ii,1} = BOS.lat.L{ii,1} - COMpos(ind1:ind2,2); % y is + to the left so COM will be stable if more to the right
        VCOM.lat.L{ii,1} = COMvel(ind1:ind2,2);

        COM.ant.L{ii,1} = BOS.ant.L{ii,1} - COMpos(ind1:ind2,1);

        % AP XCOM will differ with perturbations (depends on TM speed)
        if perts == 0 % no perturbations so use normal treadmill speeds regardless
            XCOM.ant.L{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.SSWS)./w; % relative to TM
            VCOM.ant.L{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.SSWS;
        elseif perts == 1 % perturbations so speed is affected
            % check if this is perturbation step
            if ismember(ii,pert_events.L)
                XCOM.ant.L{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.pert)./w; % relative to TM perturbation speed
                VCOM.ant.L{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.pert;
            else % use normal treadmill speeds
                XCOM.ant.L{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.SSWS)./w; % relative to TM
                VCOM.ant.L{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.SSWS;
            end
        end

        % Calculate MOS = BOS - XCOM where MOS is + if COM is within BOS
        % X is positive forward
        MOS.ant.L{ii,1} = BOS.ant.L{ii,1} - XCOM.ant.L{ii,1};
        % Y is positive to the left
        MOS.lat.L{ii,1} = BOS.lat.L{ii,1} - XCOM.lat.L{ii,1};
    end

    % RIGHT
    for ii = 1:length(RHS_index)
        % From HS to TO (ipsilateral)
        ind1 = RHS_index(ii);
        ind2 = RTO_index(find(RTO_index>ind1,1,'first'));

        % Assign BOS based on anterior/lateral foot positions
        BOS.ant.R{ii,1} = AntFoot.R(ind1:ind2,1); % X is AP
        BOS.lat.R{ii,1} = LatFoot.R(ind1:ind2,2); % Y is ML

        % Calculate XCOM = COM + VCOM/w
        XCOM.lat.R{ii,1} = COMpos(ind1:ind2,2) + COMvel(ind1:ind2,2)./w;
        COM.lat.R{ii,1} = COMpos(ind1:ind2,2) - BOS.lat.R{ii,1}; % y is + to the L so COM will be + and stable if > BOS position
        VCOM.lat.R{ii,1} = COMvel(ind1:ind2,2);
        
        COM.ant.R{ii,1} = BOS.ant.R{ii,1} - COMpos(ind1:ind2,1); % BOS should be in front of COM

        % AP XCOM will differ with perturbations (depends on TM speed)
        if perts == 0 % no perturbations so use normal treadmill speeds regardless
            XCOM.ant.R{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.SSWS)./w; % relative to TM
            VCOM.ant.R{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.SSWS;
        elseif perts == 1 % perturbations so speeed is affected
            % check if this is a perturbation step
            if ismember(ii,pert_events.R)
                XCOM.ant.R{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.pert)./w; % relative to TM
                VCOM.ant.R{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.pert;
            else
                XCOM.ant.R{ii,1} = COMpos(ind1:ind2,1) + (COMvel(ind1:ind2,1)+TMspeeds.SSWS)./w; % relative to TM
                VCOM.ant.R{ii,1} = COMvel(ind1:ind2,1)+TMspeeds.SSWS;
            end
        end

        % Calculate MOS = BOS - XCOM where MOS is + if COM is within BOS
        % X is positive forward
        MOS.ant.R{ii,1} = BOS.ant.R{ii,1} - XCOM.ant.R{ii,1};
        % Y is positive to the left
        MOS.lat.R{ii,1} = XCOM.lat.R{ii,1} - BOS.lat.R{ii,1};
    end

    % Get minimum MOS, MOS at HS, and MOS resampled per step
    for ii = 1:length(MOS.ant.L)
        MOS.ant.min.L(ii,1) = min(MOS.ant.L{ii});
        MOS.lat.min.L(ii,1) = min(MOS.lat.L{ii});
        MOS.ant.HS.L(ii,1) = MOS.ant.L{ii}(1);
        MOS.lat.HS.L(ii,1) = MOS.lat.L{ii}(1);
        MOS.ant.resample.L(1:100,ii) = resample(MOS.ant.L{ii},100,length(MOS.ant.L{ii}));
        MOS.lat.resample.L(1:100,ii) = resample(MOS.lat.L{ii},100,length(MOS.lat.L{ii}));

        COM.ant.HS.L(ii,1) = COM.ant.L{ii}(1);
        COM.lat.HS.L(ii,1) = COM.lat.L{ii}(1);

        VCOM.ant.HS.L(ii,1) = VCOM.ant.L{ii}(1);
        VCOM.lat.HS.L(ii,1) = VCOM.lat.L{ii}(1);

        XCOM.ant.HS.L(ii,1) = XCOM.ant.L{ii}(1);
        XCOM.lat.HS.L(ii,1) = XCOM.lat.L{ii}(1);
    end
    for ii = 1:length(MOS.ant.R)
        MOS.ant.min.R(ii,1) = min(MOS.ant.R{ii});
        MOS.lat.min.R(ii,1) = min(MOS.lat.R{ii});
        MOS.ant.HS.R(ii,1) = MOS.ant.R{ii}(1);
        MOS.lat.HS.R(ii,1) = MOS.lat.R{ii}(1);
        MOS.ant.resample.R(1:100,ii) = resample(MOS.ant.R{ii},100,length(MOS.ant.R{ii}));
        MOS.lat.resample.R(1:100,ii) = resample(MOS.lat.R{ii},100,length(MOS.lat.R{ii}));

        COM.ant.HS.R(ii,1) = COM.ant.R{ii}(1);
        COM.lat.HS.R(ii,1) = COM.lat.R{ii}(1);

        VCOM.ant.HS.R(ii,1) = VCOM.ant.R{ii}(1);
        VCOM.lat.HS.R(ii,1) = VCOM.lat.R{ii}(1);

        XCOM.ant.HS.R(ii,1) = XCOM.ant.R{ii}(1);
        XCOM.lat.HS.R(ii,1) = XCOM.lat.R{ii}(1);
    end

    % Sort temporally
    [~,sort_ind] = sort([LHS_index; RHS_index]);

    temp_MOS_ant_min = [MOS.ant.min.L; MOS.ant.min.R];
    MOS.sorted.ant.min(:,1) = temp_MOS_ant_min(sort_ind);

    temp_MOS_ant_HS = [MOS.ant.HS.L; MOS.ant.HS.R];
    MOS.sorted.ant.HS(:,1) = temp_MOS_ant_HS(sort_ind);
    
    temp_MOS_lat_min = [MOS.lat.min.L; MOS.lat.min.R];
    MOS.sorted.lat.min(:,1) = temp_MOS_lat_min(sort_ind);
    
    temp_MOS_lat_HS = [MOS.lat.HS.L; MOS.lat.HS.R];
    MOS.sorted.lat.HS(:,1) = temp_MOS_lat_HS(sort_ind);

    temp_COM_ant_HS = [COM.ant.HS.L; COM.ant.HS.R];
    COM.sorted.ant.HS(:,1) = temp_COM_ant_HS(sort_ind);

    temp_COM_lat_HS = [COM.lat.HS.L; COM.lat.HS.R];
    COM.sorted.lat.HS(:,1) = temp_COM_lat_HS(sort_ind);

    temp_VCOM_ant_HS = [VCOM.ant.HS.L; VCOM.ant.HS.R];
    VCOM.sorted.ant.HS(:,1) = temp_VCOM_ant_HS(sort_ind);

    temp_VCOM_lat_HS = [VCOM.lat.HS.L; VCOM.lat.HS.R];
    VCOM.sorted.lat.HS(:,1) = temp_VCOM_lat_HS(sort_ind);

    temp_XCOM_ant_HS = [XCOM.ant.HS.L; XCOM.ant.HS.R];
    XCOM.sorted.ant.HS(:,1) = temp_XCOM_ant_HS(sort_ind);

    temp_XCOM_lat_HS = [XCOM.lat.HS.L; XCOM.lat.HS.R];
    XCOM.sorted.lat.HS(:,1) = temp_XCOM_lat_HS(sort_ind);

end