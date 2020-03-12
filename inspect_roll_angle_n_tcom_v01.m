close all; clearvars; clc

%% load aligned data
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

%%
nsub = length(All_selected_sub);
All_info_trial = cell(nsub, 2);
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    beh = load(fullfile(All_path, All_filelist(All_i).name));
    nep = length(beh.file_list);
    info_trial = cell(nep, 5);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
        info_trial{ep, 3} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        info_trial{ep, 4} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
    end
    
    info_trial(:, 5) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 6) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    All_info_trial{All_i, 1} = sub_id;
    All_info_trial{All_i, 2} = info_trial;
end

%% plot for each sub
sub = 1;
info_trial = All_info_trial{sub, 2};
n_trial = length(info_trial);
% get block index
[~, ~, tmp_block_id] = unique({info_trial{:, 2}});
tmp_jump = diff(tmp_block_id);
tmp_b = 1;
b = 1;
for ep = 1:n_trial - 1
    if abs(tmp_jump(ep)) > 0
        tmp_block{b, 1} = tmp_b:ep;
        tmp_b = ep + 1;
        b = b + 1;
    end
end

f = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
subplot 211
hold on;
for ep = 1:length(tmp_block)
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = '-ro';
        case 'TR'
            linespec = '-bx';
        case 'PT'
            linespec = '-ko';
    end
    plot([info_trial{tmp_block{ep, 1}, 1}], [info_trial{tmp_block{ep, 1}, 6}], linespec);
end
xlim([0, n_trial + 1])
ylim([-600, 600])
hline([-395, 395], 'r--')
hold off
ylabel('Tcom (Nmm)');
title(All_info_trial{sub, 1});
set(gca, 'FontSize', 24);

subplot 212
hold on;
for ep = 1:length(tmp_block)
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = '-ro';
        case 'TR'
            linespec = '-bx';
        case 'PT'
            linespec = '-ko';
    end
    plot([info_trial{tmp_block{ep, 1}, 1}], abs([info_trial{tmp_block{ep, 1}, 5}]), linespec);
end
xlim([0, n_trial + 1])
ylim([0, 18])

ylabel('peak roll ({\circ})');
xlabel('trial ID');
set(gca, 'FontSize', 24);