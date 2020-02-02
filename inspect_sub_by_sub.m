close all; clearvars; clc

%% load aligned data
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/behavior/results/SXXX_temp_result.mat')
All_path = uigetdir;
All_filelist = dir(fullfile(All_path, '*_temp_result.mat'));

disp([num2cell((1:length(All_filelist))'), {All_filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');

if isempty(All_selected_sub)
    All_selected_sub = 1:length(All_filelist);
end

%%
All_dt = 0.001 * 4; % in milisecond
% color theme
All_colorset = gray;%flipud(gray);%parula;%cool;%varycolor(19);
All_ntrial = 19;
All_color_id = ceil( (1/All_ntrial) * 0.75 * length(All_colorset) );

All_cond = {'IL', 'TR', 'PT1'};

for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(1:4);
    disp(['Start processing ', sub_id, ' ...']);
    
    load(fullfile(All_path, All_filelist(All_i).name));
    
    tmp_tid_IL = [];
    tmp_tid_TR = [];
    tmp_tid_PT1 = [];
    for i = 1:length(file_list)
        if strcmp(file_list(i, 1).name(11), 'I')
            tmp_tid_IL = [tmp_tid_IL; i];
        elseif strcmp(file_list(i, 1).name(11), 'T')
            tmp_tid_TR = [tmp_tid_TR; i];
        elseif (strcmp(file_list(i, 1).name(11), 'P') && str2double(file_list(i, 1).name(16:18)) == 1)
            tmp_tid_PT1 = [tmp_tid_PT1; i];
        end
    end
    tmp_cond = {tmp_tid_IL, tmp_tid_TR, tmp_tid_PT1};
    
    ncond = length(tmp_cond);
    figure
    for i = 1:ncond
        tmp_ntrial = length(tmp_cond{i});
        for j = 1:tmp_ntrial
            tmp_onset_id = ind_lft_onset{tmp_cond{i}(j), 'h3_mm'};
            tmp_height = obj_height{tmp_cond{i}(j), 1};
            tmp_roll = angTilt2R{tmp_cond{i}(j), 1};
            
            % align all trials at the lift-onset
            tmp_trial_length = length(tmp_height);
            time_label = ( -tmp_onset_id + 1:(tmp_trial_length - tmp_onset_id) )'* All_dt;
            
            subplot(3, 2, i * 2 - 1)
            hold on
            plot(time_label, tmp_height, 'Color', All_colorset(j * All_color_id, :))
            if j == tmp_ntrial
                xlim([-2, 7])
                ylim([-100, 350])
                vline(0, '--r', 'onset')
                hline(3, ':b')
                ylabel({All_cond{i}, 'obj height (mm)'})
                if i ~= ncond
                    xticklabels([])
                else
                    xlabel('time (s)')
                end
                set(gca, 'FontSize', 18)
            end
            
            subplot(3, 2, i * 2)
            hold on
            plot(time_label, tmp_roll, 'Color', All_colorset(j * All_color_id, :))
            if j == tmp_ntrial
                xlim([-2, 7])
                vline(0, '--r', 'onset')
                ylabel('obj roll ({\circ})')
                if i ~= ncond
                    xticklabels([])
                else
                    xlabel('time (s)')
                end
                set(gca, 'FontSize', 18)
            end
        end
    end
    mtit(['Sub ', file_list(1, 1).name(2:4)], 'FontSize', 24)
    savefig(fullfile(All_path, [sub_id, '_peakRolls.fig']));
end
