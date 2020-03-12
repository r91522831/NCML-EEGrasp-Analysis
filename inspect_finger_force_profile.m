close all; clearvars; clc

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
All_info_trial = cell(nsub, 2);
All_sub_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(end-1:end);
    disp(['Start processing sub-', sub_id, ' ...']);
    
    beh = load(fullfile(All_path, All_filelist(All_i).name, 'beh', ['S0', sub_id, '_temp_result.mat']));
    nep = length(beh.file_list);
    info_trial = cell(nep, 7);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
        info_trial{ep, 3} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        info_trial{ep, 4} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
        f_mag_th = sqrt(beh.finger_Th{ep, 1}.fy .^2 + beh.finger_Th{ep, 1}.fz .^2); % F_mag_TH
        f_ang_th = atan2d(beh.finger_Th{ep, 1}.fy, beh.finger_Th{ep, 1}.fz); % F_ang_TH
        f_mag_vf = sqrt(beh.finger_V{ep, 1}.fy .^2 + beh.finger_V{ep, 1}.fz .^2); % F_mag_VF
        f_ang_vf = atan2d(beh.finger_V{ep, 1}.fy, beh.finger_V{ep, 1}.fz); % F_ang_VF
        
        d_copy_thvf = beh.finger_Th{ep, 1}.COPy - beh.finger_V{ep, 1}.COPy; % delta COPy TH - VF
        
        info_trial{ep, 5} = [f_mag_th, f_ang_th, f_mag_vf, f_ang_vf, d_copy_thvf];
    end
    
    info_trial(:, 6) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 7) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    All_info_trial{All_sub_i, 1} = ['sub-', sub_id];
    All_info_trial{All_sub_i, 2} = info_trial;
    All_sub_i = All_sub_i + 1;
end

%% plot force magnitude, force angle, and delta COPy for TH and VF in IL, TR, PT conditions

