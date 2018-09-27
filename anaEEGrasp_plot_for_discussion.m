close all; clearvars; clc

%% load aligned data
[filename, pathname, ~] = uigetfile;
load(fullfile(pathname, filename));

%% lowpass filter all data
% % % dt = diff(data{1}{1:2, 1}) * 0.001; % in second
% % % cutoff = 30; % in Hz
% % % data_filtered = data;
% % % for i = 1:length(file_list)
% % %     data_filtered{i}{:, 3:end} = filtmat_class( dt, cutoff, data{i}{:, 3:end} );
% % % end

%%
timestamp_go2fini = cell(1, 1);
obj_go2fini = cell(1, 1);
th_go2fini = cell(1, 1);
re_go2fini = cell(1, 1);
tmp_ind = 1;
for i = 1:length(file_list)
    % get only trials from the initial learing session
    if file_list(i).name(:, 11) ~= 'I'
        continue;
    end
    % plot obj
    onset_t = info_time_trigger{i, 1}(ind_lft_onset(i, 2), 1); % in ms
    go2fini = ismember(info_time_trigger{i, 2}, 2:4);
    timestamp_go2fini{tmp_ind, 1} = (info_time_trigger{i, 1}(go2fini, 1) - onset_t) * 0.001; % in s
    obj_go2fini{tmp_ind, 1} = obj{i, 1}(go2fini, :);
    th_go2fini{tmp_ind, 1} = fgr_all{i, 1}.Th{:}(go2fini, :);
    re_go2fini{tmp_ind, 1} = fgr_all{i, 1}.Re{:}(go2fini, :);
    
    tmp_ind = tmp_ind + 1;
end

% for plot
colorset = flipud(gray);%parula;%cool;%varycolor(19);
ntrial = length(timestamp_go2fini);
color_id = floor( (1/ntrial) * 0.75 * length(colorset) );

%% obj in IL
figure(1);
subplot 321
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, obj_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('x (mm)')
xticklabels([])
subplot 323
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, obj_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('y (mm)')
xticklabels([])
subplot 325
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, obj_go2fini{i}{:, 'z'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('z (mm)')
xlabel('time (s)')
subplot 322
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, rad2deg(obj_go2fini{i}{:, 'pitch'}), 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('pitch ({\circ})')
xticklabels([])
subplot 324
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, rad2deg(obj_go2fini{i}{:, 'yaw'}), 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('yaw ({\circ})')
xticklabels([])
subplot 326
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, rad2deg(obj_go2fini{i}{:, 'roll'}), 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('roll ({\circ})')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' object pose in IL'])
savefig(fullfile(pathname, 'plots', [file_list(1, 1).name(1:4), '_obj_pose_il']))

%% thumb in IL
figure(2);
subplot 331
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fx (N)')
xticklabels([])
subplot 334
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fy'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fy (N)')
xticklabels([])
subplot 337
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fz (N)')
xlabel('time (s)')
subplot 332
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'mx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('mx (N-mm)')
xticklabels([])
subplot 335
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'my'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('my (N-mm)')
xticklabels([])
subplot 338
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'mz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('mz (N-mm)')
xlabel('time (s)')
subplot 333
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('COPx (mm)')
xticklabels([])
subplot 336
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('COPy (mm)')
xticklabels([])
subplot 339
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'z'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('COPz (mm)')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' thumb in IL'])
savefig(fullfile(pathname, 'plots', [file_list(1, 1).name(1:4), '_th_il']))

%% resultant in IL
figure(3);
subplot 321
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fx (N)')
xticklabels([])
subplot 323
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fy'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fy (N)')
xticklabels([])
subplot 325
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('fz (N)')
xlabel('time (s)')
subplot 322
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'mx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('mx (N-mm)')
xticklabels([])
subplot 324
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'my'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('my (N-mm)')
xticklabels([])
subplot 326
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'mz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
ylabel('mz (N-mm)')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' resultant in IL'])
savefig(fullfile(pathname, 'plots', [file_list(1, 1).name(1:4), '_re_il']))

%% obj in TR

%% obj in PT



