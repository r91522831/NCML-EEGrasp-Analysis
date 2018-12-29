close all; clear; clc;

% combine all data
base_folder = uigetdir('Select base folder');
% folder_list = dir(fullfile(base_folder, 'S*'));

%% z_sub and g_sub into ALL_z and ALL_g
z_sub_list = dir(fullfile(base_folder, '*_z_sub.mat'));
g_sub_list = dir(fullfile(base_folder, '*_g_sub.mat'));

nb_sub = length(z_sub_list);

ALL_z = [];
ALL_g = [];

tmp_z = cell(nb_sub, 1);
tmp_g = cell(nb_sub, 1);
for sub_i = 1:nb_sub
    subID = z_sub_list(sub_i).name(1:4);
    clear g_sub;
    load(fullfile(base_folder, g_sub_list(sub_i).name));
    tmp_all_g = zeros(length(g_sub), nb_sub);
    tmp_all_g(:, sub_i) = 1;
    ALL_g = [ALL_g; tmp_all_g, g_sub];
    
    clear z_sub;
    load(fullfile(base_folder, z_sub_list(sub_i).name));
    ALL_z = [ALL_z; z_sub];
end

save(fullfile(base_folder, 'ALL_z'), 'ALL_z')
save(fullfile(base_folder, 'ALL_g'), 'ALL_g')

%% h_sub
clear h_sub;
load(fullfile(base_folder, 'h_sub.mat'));
ALL_h = h_sub;
save(fullfile(base_folder, 'ALL_h'), 'ALL_h')

%%


