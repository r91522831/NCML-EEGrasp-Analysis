close all; clear; clc;

% combine all data
base_folder = uigetdir('Select base folder');
% folder_list = dir(fullfile(base_folder, 'S*'));

%% z_sub and g_sub into ALL_z and ALL_g
z_sub_list = dir(fullfile(base_folder, '*_z_sub.mat'));
g_sub_list = dir(fullfile(base_folder, '*_g_sub.mat'));

nb_sub = length(z_sub_list);

z_cell = cell(nb_sub, 1);
g_cell = cell(nb_sub, nb_sub);

tmp_g = cell(nb_sub, 1);
for sub_i = 1:nb_sub
    subID = z_sub_list(sub_i).name(1:4);
    clear g_sub;
    load(fullfile(base_folder, g_sub_list(sub_i).name));
    
    g_cell{sub_i, sub_i} = g_sub;
    
    clear z_sub;
    load(fullfile(base_folder, z_sub_list(sub_i).name));
    z_cell{sub_i, 1} = z_sub;
end

for i = 1:sub_i
    [n, p] = size(g_cell{i, i});
    for j = 1:sub_i
        if i ~= j
            g_cell{i, j} = zeros(n, p);
        end
    end
end

ALL_g = cell2mat(g_cell);
ALL_z = cell2mat(z_cell);

save(fullfile(base_folder, 'ALL_z'), 'ALL_z')
save(fullfile(base_folder, 'ALL_g'), 'ALL_g')

%% h_sub
% h_sub_list = dir(fullfile(base_folder, '*_h_sub.mat'));
%{
bin_chan = 63;
bin_time = 40;
bin_freq = 7;
h_cell = cell(bin_chan, bin_chan);
clear h_sub;
load(fullfile(base_folder, h_sub_list(1).name)); % h_sub are the same for all subjects
tmp = h_sub(1:(bin_freq * bin_time), (bin_chan + 1):(bin_chan + bin_freq));
for i = 1:bin_chan
    h_cell{i, i} = tmp;
end

for i = 1:bin_chan
    [m, q] = size(h_cell{i, i});
    for j = 1:bin_chan
        if i ~= j
            h_cell{i, j} = zeros(m, q);
        end
    end
end
ALL_h = cell2mat(h_cell);

% % % clear h_sub;
% % % load(fullfile(base_folder, 'h_sub.mat'));
% % % ALL_h = h_sub;

save(fullfile(base_folder, 'ALL_h'), 'ALL_h')
%}


%%


