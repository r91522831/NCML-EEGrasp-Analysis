close all; clearvars; clc

%% load aligned data

All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
%%
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end
    
% load(fullfile(All_path, 'behavior_pRoll_23subs.mat'));
load(fullfile(All_path, 'behavior_pRoll_12subs.mat'));

All_period_3mm2objweight = nan(length(All_filelist), 1);
for All_i = All_selected_sub
    clearvars -except All_* all_sub_mx_pRoll; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    
    All_period_3mm2objweight(All_i) = nanmean(abs(ind_lft_onset{:, 3} - ind_lft_onset{:, 6}));
    
    
end