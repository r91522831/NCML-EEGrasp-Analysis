close all; clearvars; clc

%% load aligned data

All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
%%
All_selected_sub = input('Which subject(s) to process? ');

% [filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/matlab data/sandbox/SXXX_fingerInfo.mat')

load(fullfile(All_path, 'behavior_pRoll_23subs.mat'));

for All_i = All_selected_sub
    clearvars -except All_* all_sub_mx_pRoll; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    
    %%
    plot_trial = input('Which trial(s) to plot? ');
    if isempty(plot_trial)
        plot_trial = 1:length(resultantF);
    end
    
    for t = plot_trial%1:length(resultantF)
        non_zero = [true(1); info_time_trigger{t, 1}(2:end) ~= 0];
        figure
        for c = 1:6
            subplot(3, 3, c)
            plot(info_time_trigger{t, 1}(non_zero)./1000, resultantF{t, 1}{non_zero, c});
            vline(info_time_trigger{t, 1}(ind_lft_onset{t, 3})./1000);
        end
        subplot(3, 3, 7)
        plot(info_time_trigger{t, 1}(non_zero)./1000, obj_height{t, 1}(non_zero));
        vline(info_time_trigger{t, 1}(ind_lft_onset{t, 3})./1000);
        subplot(3, 3, 8)
        tmp_ang = angTilt2R{t, 1};
        for j = 1:length(angTilt2R{t, 1})
            if angTilt2R{t, 1}(j) > 90
                tmp_ang(j) = angTilt2R{t, 1}(j) - 180;
            elseif angTilt2R{t, 1}(j) < -90
                tmp_ang(j) = angTilt2R{t, 1}(j) + 180;
            end
        end
        plot(info_time_trigger{t, 1}(non_zero)./1000, tmp_ang(non_zero));
        vline(info_time_trigger{t, 1}(ind_lft_onset{t, 3})./1000);

        suptitle(['Trial ', num2str(t)])
        set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);
    end

end