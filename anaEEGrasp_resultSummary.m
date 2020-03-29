close all; 
clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, 'sub*'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

%%
nsub = length(All_selected_sub);
All_info_trial = cell(nsub, 3);
All_sub_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(end-1:end);
    disp(['Start processing sub-', sub_id, ' ...']);
    
    beh = load(fullfile(All_path, All_filelist(All_i).name, 'beh', 'mat', ['S0', sub_id, '_temp_result.mat']));
    
    dt = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % dt in ms
    tmp_win = round(-1000/dt):round(2000/dt); % 1 s before lift onset to 2 s after lift onset
    
    nep = length(beh.file_list);
    info_trial = cell(nep, 14);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
        info_trial{ep, 3} = beh.file_list(ep).name(20); % handle side
        if length(beh.file_list(ep).name) >= 25
            info_trial{ep, 4} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        end
        info_trial{ep, 5} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
        deltaCOPy = beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'COPy'} - beh.finger_V{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'COPy'};
        magFth = sqrt(beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'}.^2 + beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fz'}.^2);
        deltaFy = beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'} - beh.finger_V{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'};
        gripF = 0.5 * abs(beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fz'} - beh.finger_V{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fz'});
        
        info_trial{ep,  8} = deltaCOPy;
        info_trial{ep,  9} = magFth;
        info_trial{ep, 10} = deltaFy;
        info_trial{ep, 11} = gripF;
        
        f_mag_th = sqrt(beh.finger_Th{ep, 1}.fy .^2 + beh.finger_Th{ep, 1}.fz .^2); % F_mag_TH
        f_ang_th = atan2(beh.finger_Th{ep, 1}.fy, beh.finger_Th{ep, 1}.fz); % F_ang_TH
        f_mag_vf = sqrt(beh.finger_V{ep, 1}.fy .^2 + beh.finger_V{ep, 1}.fz .^2); % F_mag_VF
        f_ang_vf = atan2(beh.finger_V{ep, 1}.fy, beh.finger_V{ep, 1}.fz); % F_ang_VF
        
        f_ang_vf(f_ang_vf < 0) = f_ang_vf(f_ang_vf < 0) + 2 * pi; % prevent jumping of angle close to 180 degree

        d_copy_thvf = beh.finger_Th{ep, 1}.COPy - beh.finger_V{ep, 1}.COPy; % delta COPy TH - VF
        d_fy_thvf = beh.finger_Th{ep, 1}.fy - beh.finger_V{ep, 1}.fy;
        grip_f = 0.5 * ( abs(beh.finger_Th{ep, 1}.fz) + abs(beh.finger_V{ep, 1}.fz) );
        f_n_th = beh.finger_Th{ep, 1}.fz;
        f_n_vf = beh.finger_V{ep, 1}.fz;
        f_y_th = beh.finger_Th{ep, 1}.fy;
        f_y_vf = beh.finger_V{ep, 1}.fy;
        
        profile = array2table( [f_mag_th, f_ang_th, f_mag_vf, f_ang_vf, d_copy_thvf, d_fy_thvf, f_n_th, f_n_vf, f_y_th, f_y_vf, grip_f], ...
                                          'VariableNames', {'magFth', 'angFth', 'magFvf', 'angFvf', 'dCOPy', 'dFy', 'FnTH', 'FnVF', 'FyTH', 'FyVF', 'GripF'} );
        info_trial{ep, 12} = profile;
        time_all = beh.info_onset_time(ep, :);
        info_trial{ep, 13} = time_all;
        
        lft_ind = time_all{1, 'lft_ind'};
        lft_time = time_all{1, 'lft_time'};
        lft_win = ( lft_ind + tmp_win );
        lft_win_time = tmp_win * dt * 0.001; % in s
        lft_win_time = array2table(lft_win_time', 'VariableNames', {'time'});
        info_trial{ep, 14} = horzcat(lft_win_time, profile(lft_win, :));
    end
    
    info_trial(:, 6) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 7) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    info_trial = cell2table(info_trial, 'VariableNames', {'trial', 'context', 'side', 'extraW', 'inBlock', 'pRoll', 'Mcom', 'dCOPy', 'magFth', 'dFy', 'gripF', 'profile', 'lftOnset', 'onset_aligned'});
    
    All_info_trial{All_sub_i, 1} = ['sub-', sub_id];
    All_info_trial{All_sub_i, 2} = info_trial;
    All_info_trial{All_sub_i, 3} = dt;
    All_sub_i = All_sub_i + 1;
end
All_info_trial = cell2table(All_info_trial, 'VariableNames', {'subID', 'data', 'dt'});
save(fullfile(All_path, 'behavior_resultSummary.mat'), 'All_info_trial');