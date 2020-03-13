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
All_peakRoll = nan(95, length(All_selected_sub));
All_TCOM = nan(95, 2, length(All_selected_sub));
All_ind = 1;
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name;
    disp(['Start processing ', sub_id, ' ...']);
    
    sub_n = sub_id(end-1:end);
    filename = fullfile(All_path, All_filelist(All_i).name, 'beh', 'mat', ['S0', sub_n, '_temp_result.mat']);
    beh = load(filename);
    
    trial_id = cell2mat(cellfun(@(x) str2double(x(7:9)), {beh.file_list.name}, 'UniformOutput', false));
    for ep = 1:length(trial_id)
        All_TCOM(trial_id(ep), :, All_ind) = [beh.peak_mx{ep, 'peakMx'}, sign(beh.peak_mx{ep, 'peakMx'}) * 395];
        All_peakRoll(trial_id(ep), All_ind) = abs(beh.peak_roll{ep, 'peakRoll'});
    end
    All_ind = All_ind + 1;
end

TCOM = All_TCOM;
peakRoll = All_peakRoll;