close all; clear; clc;

% combine all data
base_folder = uigetdir('Select base folder');
folder_list = dir(fullfile(base_folder, 'S*'));

%%
cond_names = {'IL', 'TR', 'PT1', 'PT2', 'PT3'};
cond_nb = length(cond_names);
% typeproc - type of processing: 1 process the raw channel data
%                                0 process the ICA component data
typeproc = str2double(input('\nChoose type of processing [1: raw, 0: component]: (default: 0) ', 's'));

ALL_tf = cell(length(folder_list), 2);
nb_sub = length(folder_list);
ALL_z = [];
ALL_g = [];
for sub_i = 1:2%nb_sub
    subID = folder_list(sub_i).name;
    
    
    
    
    ALL_z = [ALL_z; z_sub];
    tmp_all_g = zeros(EEG.trials, nb_sub);
    tmp_all_g(:, sub_i) = 1;
    ALL_g = [ALL_g; tmp_all_g, g_sub];
    ALL_h = h_sub;
end
%%
