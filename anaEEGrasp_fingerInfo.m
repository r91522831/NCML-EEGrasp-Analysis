close all; clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_aligned_data.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

% [filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/preliminary results/SXXX_aligned_data.mat')

All_cutoff = 5; % lowpass cutoff freq for behavior, in Hz

%%
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    % file_list is included in the SXXX_aligned_data.mat file.
    % file_list contains all behavior files name with wrong directory.

    n_PSonObj = 8;
    n_PSonFgr = 3;
    kinetic_channels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};

    %% 
    resultantF = cell(size(file_list));
    finger_Th = cell(size(file_list));
    finger_V = cell(size(file_list));
    info_time_trigger = cell(length(file_list), 2);
    missing_info = cell(length(file_list), 1);
    coord_obj = cell(size(file_list));
    fgr_on_obj = cell(size(file_list));
    fgr_on_obj_kinetic = cell(size(file_list));

    %% lowpass filtered
    % % % dt = 4/1000; % unit: second
    % % % for i = 1:length(mx)
    % % %     mx_filtered{i} = filtmat_class( dt, cutoff_plot, mx{i} );
    % % %     obj_height_filtered{i} = filtmat_class( dt, cutoff_plot, obj_height{i} );
    % % % end


    %--------------------------------------------------------------------------
    % assign filtered or raw data to analyze
    input = data;
    input_surface = data_aligned2surface;

    %%
    for i = 1:length(file_list)
        % get which side of the handle on the object is grasped
        obj_side = file_list(i).name(end-4);
        % get time stamp
        tmp_tstamp = input{i}{:, {'time_ms'}};
        info_time_trigger{i, 1} = tmp_tstamp;
        % get trigger indices
        tmp_audio = input{i}{:, {'trigger'}};
        info_time_trigger{i, 2} = tmp_audio;
        
        tmp_audio = info_time_trigger{i, 2};
        dt = 0.001 * diff(tmp_tstamp(1:2)); % in second

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Check the PS sensor in each frame
        % lowpass markers coordinates without touching conditions
        markers = cell(n_PSonObj, 1);
        for m = 1:n_PSonObj
            lp_coord = filtmat_class( dt, All_cutoff, input{i}{:, var_PS{m}});
            markers{m, 1} = [lp_coord, input{i}{:, var_PS_cond{m}}];
        end
        markers_fgr = cell(n_PSonFgr, 1);
        for m = 1:n_PSonFgr
            lp_coord = filtmat_class( dt, All_cutoff, input{i}{:, var_PS{m + 8}});
            markers_fgr{m, 1} = [lp_coord, input{i}{:, var_PS_cond{m + 8}}];
        end

        tmp_fgr_on_obj = cell(n_PSonFgr, 1);
        tmp_missing = 1;
        for time_id = 1:height(input{i})
            % for markers on the object
            tmp_markers = cell(n_PSonObj, 1);
            tmp_cond = zeros(n_PSonObj, 1);
            for m = 1:n_PSonObj
                tmp_markers{m, 1} = markers{m, 1}(time_id, 1:3);
                tmp_cond(m, 1) = markers{m, 1}(time_id, 4);
            end
            tmp_yz_markers = [2, 4, 5]; % default markers
            tmp_xz_markers = [3, 6, 8];
            tmp_origin_markers = 1:8;

            % for markers on three finger tips
            tmp_markers_fgr = cell(n_PSonFgr, 1);
            tmp_cond_fgr = zeros(n_PSonFgr, 1);
            for m = 1:n_PSonFgr
                tmp_markers_fgr{m, 1} = markers_fgr{m, 1}(time_id, 1:3);
                tmp_cond_fgr(m, 1) = markers_fgr{m, 1}(time_id, 4);
            end

            %% might not need to check the missing markers, since it was done when converting from .tdms to .cvs
            % this will introduce extra noise
            %{
            % if there is any marker missing in this frame
            if any(tmp_cond < 0) || any(tmp_cond_fgr < 0)
                tmp_missing_markers = find(tmp_cond < 0);
    %             tmp_missing_markers_fgr = find(tmp_cond_fgr < 0);
                % check if there is at least 3 within [1, 2, 4, 5]
                tmp_critical = [1, 2, 4, 5];
                tmp_missed = ismember(tmp_critical, tmp_missing_markers);


                switch sum(tmp_missed)
                    case {0, 1}
                        tmp_yz_markers = tmp_critical( ~tmp_missed );
                    case 2
                        tmp_ind_missed = find(tmp_missed);
                        tmp_not_picked = tmp_missed;
                        tmp_not_picked(tmp_ind_missed(1)) = 0;
                        tmp_yz_markers = tmp_critical( ~tmp_not_picked );
                        disp(['trial ', num2str(i), ' missed 2 markers for yz plan!']);
                    case 3
                        tmp_ind_missed = find(tmp_missed);
                        tmp_not_picked = tmp_missed;
                        tmp_not_picked(tmp_ind_missed(1:2)) = 0;
                        tmp_yz_markers = tmp_critical( ~tmp_not_picked );
                        disp(['trial ', num2str(i), ' missed 3 markers for yz plan!']);
                    case 4
                        % use [2, 4, 5]
                        disp(['trial ', num2str(i), ' missed all markers for yz plan!']);
                end

                % check if there is all within [3, 6, 8]
                tmp_missed = ismember(tmp_xz_markers, tmp_missing_markers);
                if sum(tmp_missed) > 0
                    disp(['trial ', num2str(i), ' missed too many markers for xz plan!']);
    %                 continue; % still use the missed marker with smoothened data
                end
                % remove missing marker from markers used to get origin
                tmp_origin_markers = tmp_origin_markers(tmp_cond >= 0);

                % ------------------------------------------------------------
                % log missing marker information
                missing_info{i, 1}{tmp_missing, 1} = time_id;
                missing_info{i, 1}{tmp_missing, 2} = find(tmp_cond < 0);
                missing_info{i, 1}{tmp_missing, 3} = (tmp_cond < 0);
                missing_info{i, 1}{tmp_missing, 4} = find(tmp_cond_fgr < 0);
                missing_info{i, 1}{tmp_missing, 5} = (tmp_cond_fgr < 0);
                tmp_missing = tmp_missing + 1;
            end
            %}



            % compute coordinate bases from yz and xz plan vectors
            tmp_coordBasis = coordBasisOnObj(tmp_markers, tmp_yz_markers, tmp_xz_markers);
            % ----------------------------------------------------------------
            % reconstruct the origins (L and R) using the object rigid dimensions
            tmp_coordOrigin = coordOriginOnObj(tmp_markers, tmp_origin_markers, tmp_coordBasis, obj_side);
            coord_obj{i}{time_id, 1} = [tmp_coordBasis, tmp_coordOrigin];

            % compute the finger location in the coordinate frame fixed on the object
            for m = 1:n_PSonFgr
                tmp_fgr = (tmp_markers_fgr{m, 1} - transpose( tmp_coordOrigin{:, 'origin'} )) * table2array(tmp_coordBasis);
                tmp_fgr_on_obj{m, 1}(time_id, :) = array2table([tmp_fgr, tmp_cond_fgr(m, 1)], 'VariableNames', {'x', 'y', 'z', 'cond'});
            end
        end
        fgr_on_obj{i} = array2table(tmp_fgr_on_obj', 'VariableNames', {'Th', 'In', 'Mi'});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% get finger kinetic data
        switch obj_side
            case 'R'
                tmp_fTh = input{i}{:, var_ATI{1}};
                tmp_fV  = input{i}{:, var_ATI{3}};

                tmp_fTh_surface = input_surface{i}{:, var_ATI{1}};
                tmp_fV_surface  = input_surface{i}{:, var_ATI{3}};
            case 'L'
                tmp_fTh = input{i}{:, var_ATI{2}};
                tmp_fV  = input{i}{:, var_ATI{4}};

                tmp_fTh_surface = input_surface{i}{:, var_ATI{2}};
                tmp_fV_surface  = input_surface{i}{:, var_ATI{4}};
            otherwise
                disp('something wrong with the file name!!')
        end

        % get figer kinetics and resultant kinetics with respect to the
        % corresponding handle center
        finger_Th{i} = array2table(tmp_fTh, 'VariableNames', kinetic_channels);
        finger_V{i} = array2table(tmp_fV, 'VariableNames', kinetic_channels);
        resultantF{i} = array2table(tmp_fTh + tmp_fV, 'VariableNames', kinetic_channels);

        % --------------------------------------------------------------------
        % find th_cop_y and v_cop_y from ATI force sensors
        % finger kinetics with respect to individual contact surfaces
        tmp_finger_Th_surface = array2table(tmp_fTh_surface, 'VariableNames', kinetic_channels);
        tmp_finger_V_surface = array2table(tmp_fV_surface, 'VariableNames', kinetic_channels);
        % ATI Nano 25 with SI-125-3 calibration: Fz resolution is 1/16 N
        tmp_reslu = 1/16 * 10; % N
        % --------------------------------------------------------------------
        % the ATI data coordinate frame have been translated and rotated to the
        % center of each contact surface
        z_center2finger = 3 + 21.6 + 3 + 2; % in mm; center piece + ATI Nano25 + mounting + surface
        tmp_Th_cop = centerOfPressure(finger_Th{i}, -z_center2finger, tmp_reslu);
        tmp_V_cop = centerOfPressure(finger_V{i}, z_center2finger, tmp_reslu);

        % --------------------------------------------------------------------
        % cop z equals ATI surface to finger contact surface
        z_ATI2finger = 3 + 2; % in mm; mounting + surface
        tmp_Th_cop_surface = centerOfPressure(tmp_finger_Th_surface, -z_ATI2finger, tmp_reslu);
        tmp_V_cop_surface = centerOfPressure(tmp_finger_V_surface, z_ATI2finger, tmp_reslu);

        tmp_avgTh_cop = array2table(mean(cat(3, table2array(tmp_Th_cop), table2array(tmp_Th_cop_surface)), 3), 'VariableNames', {'COPx', 'COPy', 'COPz'});
        tmp_avgV_cop = array2table(mean(cat(3, table2array(tmp_V_cop), table2array(tmp_V_cop_surface)), 3), 'VariableNames', {'COPx', 'COPy', 'COPz'});

        finger_Th{i} = [finger_Th{i}, tmp_avgTh_cop];
        finger_V{i}  = [finger_V{i} , tmp_avgV_cop];

        tmp_fgr_on_obj_kinetic = { tmp_avgTh_cop, tmp_avgV_cop };
        fgr_on_obj_kinetic{i} = array2table(tmp_fgr_on_obj_kinetic, 'VariableNames', {'Th', 'Vi'});

        %%
        disp(i);
    end

    save(fullfile(All_path, [sub_id, '_fingerInfo.mat']), 'file_list', 'data', 'var_PS', 'var_PS_cond', 'info_time_trigger', 'missing_info', 'coord_obj', 'fgr_on_obj', 'finger_Th', 'finger_V', 'resultantF', 'fgr_on_obj_kinetic', '-v7');
    disp(['Finish processing ', sub_id, '!']);
end