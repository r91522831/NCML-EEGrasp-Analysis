close all; clearvars; clc

%% load aligned data
[filename, pathname_plot, ~] = uigetfile;
load(fullfile(pathname_plot, filename));

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
in_go2fini = cell(1, 1);
mi_go2fini = cell(1, 1);
tmp_ind = 1;
session = 'PT'; %'TR'; %'IL';
for i = 1:length(file_list)
    % get only trials from the initial learing session
    if ~strcmp(file_list(i).name(:, 11:12), session)
        continue;
    end
    
    obj_side = file_list(i).name(:, end - 4);
    % plot obj
    onset_t = info_time_trigger{i, 1}(ind_lft_onset(i, 2), 1); % in ms
    go2fini = ismember(info_time_trigger{i, 2}, 2:4);
    timestamp_go2fini{tmp_ind, 1} = (info_time_trigger{i, 1}(go2fini, 1) - onset_t) * 0.001; % in s
    obj_go2fini{tmp_ind, 1} = obj{i, 1}(go2fini, :);
    th_go2fini{tmp_ind, 1} = fgr_all{i, 1}.Th{:}(go2fini, :);
    re_go2fini{tmp_ind, 1} = fgr_all{i, 1}.Re{:}(go2fini, :);
    in_go2fini{tmp_ind, 1} = fgr_all{i, 1}.In{:}(go2fini, :);
    mi_go2fini{tmp_ind, 1} = fgr_all{i, 1}.Mi{:}(go2fini, :);
    
    tmp_ind = tmp_ind + 1;
end

ini_id = [];
if strcmp(session, 'PT')
    ini_id = 3;
    timestamp_go2fini = timestamp_go2fini(ini_id:3:end, 1);
    obj_go2fini = obj_go2fini(ini_id:3:end, 1);
    th_go2fini = th_go2fini(ini_id:3:end, 1);
    re_go2fini = re_go2fini(ini_id:3:end, 1);
    in_go2fini = in_go2fini(ini_id:3:end, 1);
    mi_go2fini = mi_go2fini(ini_id:3:end, 1);
end


% for plot
colorset = flipud(gray);%parula;%cool;%varycolor(19);
ntrial = length(timestamp_go2fini);
color_id = ceil( (1/ntrial) * 0.75 * length(colorset) );

tmp_t = zeros(ntrial, 2);
for i = 1:ntrial
    tmp_t(i, 1) = min(timestamp_go2fini{i});
    tmp_t(i, 2) = max(timestamp_go2fini{i});
end
min_t = min(tmp_t(:, 1)) - 0.05;
max_t = max(tmp_t(:, 2)) + 0.05;



%% obj in tag
figure(1);
subplot 321
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, obj_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('x (mm)')
xticklabels([])
subplot 323
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, obj_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
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
xlim([min_t, max_t])
ylabel('pitch ({\circ})')
xticklabels([])
subplot 324
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, rad2deg(obj_go2fini{i}{:, 'yaw'}), 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('yaw ({\circ})')
xticklabels([])
subplot 326
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, rad2deg(obj_go2fini{i}{:, 'roll'}), 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('roll ({\circ})')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' object pose in ', session, num2str(ini_id), ' ', obj_side])
savefig(fullfile(pathname_plot, 'plots', [file_list(1, 1).name(1:4), ['_obj_pose_', session, num2str(ini_id)]]))

%% thumb in tag
figure(2);
subplot 321
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fx (N)')
xticklabels([])
subplot 323
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fy'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fy (N)')
xticklabels([])
subplot 325
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'fz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fz (N)')
xlabel('time (s)')
subplot 322
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'mx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('mx (N-mm)')
xticklabels([])
subplot 324
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'my'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('my (N-mm)')
xticklabels([])
subplot 326
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'mz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('mz (N-mm)')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' thumb in ', session, num2str(ini_id), ' ', obj_side])
savefig(fullfile(pathname_plot, 'plots', [file_list(1, 1).name(1:4), ['_th_', session, num2str(ini_id)]]))

%% resultant in tag
figure(3);
subplot 321
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fx (N)')
xticklabels([])
subplot 323
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fy'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fy (N)')
xticklabels([])
subplot 325
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'fz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('fz (N)')
xlabel('time (s)')
subplot 322
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'mx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('mx (N-mm)')
xticklabels([])
subplot 324
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'my'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('my (N-mm)')
xticklabels([])
subplot 326
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, re_go2fini{i}{:, 'mz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylabel('mz (N-mm)')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' resultant in ', session, num2str(ini_id), ' ', obj_side])
savefig(fullfile(pathname_plot, 'plots', [file_list(1, 1).name(1:4), ['_re_', session, num2str(ini_id)]]))

%% fingers in tag
figure(4);
subplot 341
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'COPx'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-20, 20])
ylabel('COPx (mm) ATI')
xticklabels([])
subplot 345
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'COPy'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-60, 60])
ylabel('COPy (mm) ATI')
xticklabels([])
subplot 349
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'COPz'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-20, -10])
ylabel('COPz (mm) ATI')
xlabel('time (s)')
subplot 342
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-20, 20])
ylabel('COPx (mm)')
xticklabels([])
subplot 346
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-60, 60])
ylabel('COPy (mm)')
xticklabels([])
subplot(3, 4, 10)
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, th_go2fini{i}{:, 'z'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-50, -25])
ylabel('COPz (mm)')
xlabel('time (s)')
subplot 343
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, in_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-20, 20])
ylabel('IN COPx (mm)')
xticklabels([])
subplot 347
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, in_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-60, 60])
ylabel('IN COPy (mm)')
xticklabels([])
subplot(3, 4, 11)
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, in_go2fini{i}{:, 'z'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([25, 50])
ylabel('IN COPz (mm)')
xlabel('time (s)')
subplot 344
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, mi_go2fini{i}{:, 'x'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-20, 20])
ylabel('MI COPx (mm)')
xticklabels([])
subplot 348
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, mi_go2fini{i}{:, 'y'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([-60, 60])
ylabel('MI COPy (mm)')
xticklabels([])
subplot(3, 4, 12)
hold on
for i = 1:ntrial
    plot(timestamp_go2fini{i}, mi_go2fini{i}{:, 'z'}, 'Color', colorset(i * color_id, :))
    vline(0)
end
hold off
xlim([min_t, max_t])
ylim([25, 50])
ylabel('MI COPz (mm)')
xlabel('time (s)')
mtit([file_list(1, 1).name(1:4), ' fingers in ', session, num2str(ini_id), ' ', obj_side])
savefig(fullfile(pathname_plot, 'plots', [file_list(1, 1).name(1:4), ['_fgr_', session, num2str(ini_id)]]))



