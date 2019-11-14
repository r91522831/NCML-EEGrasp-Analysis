close all; clearvars; clc

%% load aligned data

All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    
    sub_id = All_filelist(All_i).name(1:4);
    tmp_ind_onset = 3; % 3 mm 
    
    fy_onset = nan(length(resultantF), 1);
    for i = 1:length(resultantF)
        fy_onset(i) = resultantF{i, 1}{ind_lft_onset{i, tmp_ind_onset}, 'fy'};
    end
    
    save(fullfile(All_path, [sub_id, '_temp_result.mat']), 'fy_onset', '-append')
end