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
        % check if the handle was rotated
        tmp_ntrial = str2double(file_list(i).name(7:9));
        if tmp_ntrial > 15 && rem(tmp_ntrial - 15, 5) == 1
            handle_rotated = true;
        else
            handle_rotated = false;
        end
        
        % get time stamp
        tmp_tstamp = input{i}{:, {'time_ms'}};
        info_time_trigger{i, 1} = tmp_tstamp;
        % get trigger indices
        tmp_audio = input{i}{:, {'trigger'}};
        info_time_trigger{i, 2} = tmp_audio;
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
            
            
            % compute coordinate bases
            if handle_rotated
                y_vector = tmp_markers{8, 1} - tmp_markers{6, 1};
                z_vector = tmp_markers{7, 1} - tmp_markers{5, 1};
                tmp_origin_markers = [6, 8, 4, 3];
            else
                y_vector = tmp_markers{2, 1} - tmp_markers{1, 1};
                z_vector = tmp_markers{5, 1} - tmp_markers{7, 1};
                tmp_origin_markers = [1, 2, 3, 4];
            end
            tmp_coord_obj_y = y_vector / norm(y_vector);
            coord_obj_z = z_vector / norm(z_vector);
            tmp = cross(tmp_coord_obj_y, coord_obj_z);
            coord_obj_x = tmp / norm(tmp);
            tmp = cross(coord_obj_z, coord_obj_x);
            coord_obj_y = tmp / norm(tmp);
            
            tmp = [coord_obj_x; coord_obj_y; coord_obj_z]';
            tmp_coordBasis = array2table(tmp, 'RowNames', {'x', 'y', 'z'}, 'VariableNames', {'x_axis', 'y_axis', 'z_axis'});
            
            % ----------------------------------------------------------------
            % reconstruct the origins (L and R) using the object rigid dimensions
            tmp_coordOrigin = TShape_coordOriginOnObj(tmp_markers, tmp_origin_markers, tmp_coordBasis);
            coord_obj{i}{time_id, 1} = [tmp_coordBasis, tmp_coordOrigin];
            
            % for markers on three finger tips
            tmp_markers_fgr = cell(n_PSonFgr, 1);
            tmp_cond_fgr = zeros(n_PSonFgr, 1);
            for m = 1:n_PSonFgr
                tmp_markers_fgr{m, 1} = markers_fgr{m, 1}(time_id, 1:3);
                tmp_cond_fgr(m, 1) = markers_fgr{m, 1}(time_id, 4);
            end
            % compute the finger location in the coordinate frame fixed on the object
            for m = 1:n_PSonFgr
                tmp_fgr = (tmp_markers_fgr{m, 1} - transpose( tmp_coordOrigin{:, 'origin'} )) * table2array(tmp_coordBasis);
                tmp_fgr_on_obj{m, 1}(time_id, :) = array2table([tmp_fgr, tmp_cond_fgr(m, 1)], 'VariableNames', {'x', 'y', 'z', 'cond'});
            end
        end
        fgr_on_obj{i} = array2table(tmp_fgr_on_obj', 'VariableNames', {'Th', 'In', 'Mi'});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% get finger kinetic data
        tmp_fTh = input{i}{:, var_ATI{1}};
        tmp_fV  = input{i}{:, var_ATI{2}}; % fx1, fy1, fz1
        
        tmp_fTh_surface = input_surface{i}{:, var_ATI{1}};
        tmp_fV_surface  = input_surface{i}{:, var_ATI{2}}; % fx1, fy1, fz1

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