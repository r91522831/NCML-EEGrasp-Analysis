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
nsub = length(All_selected_sub);
All_info_trial = cell(nsub, 3);
All_sub_i = 1;
for All_i = All_selected_sub
    clearvars -except All_*; close all;
    sub_id = All_filelist(All_i).name(end-1:end);
    disp(['Start processing sub-', sub_id, ' ...']);
    
    beh = load(fullfile(All_path, All_filelist(All_i).name, 'beh', 'mat', ['S0', sub_id, '_temp_result.mat']));
    nep = length(beh.file_list);
    info_trial = cell(nep, 8);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
% % %         info_trial{ep, 3} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        info_trial{ep, 4} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
        f_mag_th = sqrt(beh.finger_Th{ep, 1}.fy .^2 + beh.finger_Th{ep, 1}.fz .^2); % F_mag_TH
        f_ang_th = atan2(beh.finger_Th{ep, 1}.fy, beh.finger_Th{ep, 1}.fz); % F_ang_TH
        f_mag_vf = sqrt(beh.finger_V{ep, 1}.fy .^2 + beh.finger_V{ep, 1}.fz .^2); % F_mag_VF
        f_ang_vf = atan2(beh.finger_V{ep, 1}.fy, beh.finger_V{ep, 1}.fz); % F_ang_VF
        
        f_ang_vf(f_ang_vf < 0) = f_ang_vf(f_ang_vf < 0) + 2 * pi; % prevent jumping of angle close to 180 degree

        d_copy_thvf = beh.finger_Th{ep, 1}.COPy - beh.finger_V{ep, 1}.COPy; % delta COPy TH - VF
        d_fy_thvf = beh.finger_Th{ep, 1}.fy - beh.finger_V{ep, 1}.fy;
        f_n_th = beh.finger_Th{ep, 1}.fz;
        f_n_vf = beh.finger_V{ep, 1}.fz;
        f_y_th = beh.finger_Th{ep, 1}.fy;
        f_y_vf = beh.finger_V{ep, 1}.fy;
        
        info_trial{ep, 5} = [f_mag_th, f_ang_th, f_mag_vf, f_ang_vf, d_copy_thvf, d_fy_thvf, f_n_th, f_n_vf, f_y_th, f_y_vf];
        info_trial{ep, 6} = beh.info_onset_time(ep, :);
    end
    
    info_trial(:, 7) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 8) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    
    All_info_trial{All_sub_i, 1} = ['sub-', sub_id];
    All_info_trial{All_sub_i, 2} = info_trial;
    All_info_trial{All_sub_i, 3} = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % dt in ms
    All_sub_i = All_sub_i + 1;
end

%% plot force magnitude, force angle, and delta COPy for TH and VF in IL, TR, PT conditions
sub_i = 1;
sub_id = All_info_trial{sub_i, 1};
data = All_info_trial{sub_i, 2};
dt = All_info_trial{sub_i, 3}; % in ms

tmp_win = round(-1000/dt):round(2000/dt);
tch_time = nan(length(data), 1);
f_mag_th = nan(length(data), length(tmp_win));
f_mag_vf = f_mag_th; f_ang_th = f_mag_th; f_ang_vf = f_mag_th; d_copy_thvf = f_mag_th;
d_fy_thvf = f_mag_th; f_n_th = f_mag_th; f_n_vf = f_mag_th; f_y_th = f_mag_th; f_y_vf = f_mag_th;
for ep = 1:length(data)
    lft_ind = data{ep, 6}{1, 'lft_ind'};
    lft_time = data{ep, 6}{1, 'lft_time'};
    tch_time(ep, 1) = data{ep, 6}{1, 'tch_time'} - lft_time;
    lft_win = ( lft_ind + tmp_win )';
    lft_win_time = tmp_win * dt * 0.001; % in s
    
    f_mag_th(ep, :) = data{ep, 5}(lft_win, 1)';
    f_mag_vf(ep, :) = data{ep, 5}(lft_win, 3)';
    f_ang_th(ep, :) = data{ep, 5}(lft_win, 2)';
    f_ang_vf(ep, :) = data{ep, 5}(lft_win, 4)';
    d_copy_thvf(ep, :) = data{ep, 5}(lft_win, 5)';
    d_fy_thvf(ep, :) = data{ep, 5}(lft_win, 6)';
    f_n_th(ep, :) = data{ep, 5}(lft_win, 7)';
    f_n_vf(ep, :) = data{ep, 5}(lft_win, 8)';
    f_y_th(ep, :) = data{ep, 5}(lft_win, 9)';
    f_y_vf(ep, :) = data{ep, 5}(lft_win, 10)';
end

%% plot 3D plot: Fn, dFy, dCOPy
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
[~, lft_onset] = min(abs(lft_win_time));

figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
ep_plot = 1:length(tmp_block); %5:length(tmp_block);
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    x = d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset);
    y = d_fy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset);
    z = f_n_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset);
    
    plot3(x, y, z, 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 5;
        else
            fw = 'normal';
            fs = 5; %length(tmp_block{ep, 1}) - i;
        end
        x_t = d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset);
        y_t = d_fy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset);
        z_t = f_n_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset);
        text(x_t, y_t, z_t, num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
end
% % % % % % legend('IL', 'TR', 'PT');
% % % % sub-04: 18, 15, 12; sub-02: 13, 10, 7; sub-01: 18, 13, 8
text(-20, 18, 35, 'IL', 'Color', 'r', 'FontSize', 18);
text(-20, 15, 35, 'TR', 'Color', 'b', 'FontSize', 18);
text(-20, 12, 35, 'PT', 'Color', 'k', 'FontSize', 18);
grid on
set(gca, 'View', [70, 35])
axis equal
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('{\Delta}Fy_{TH-VF} (N)');
zlabel('Fn_{TH} (N)');
% % % vline(0, '--k')
% % % hline(0, '--k')
hold off

title(All_info_trial{sub, 1})











%% plot force in polar coordinate
%{
ep_ind = [strcmpi(data(:, 2), 'IL'), strcmpi(data(:, 2), 'TR'), strcmpi(data(:, 2), 'PT')];
tmp_lft_win_time = repmat(lft_win_time, length(data), 1);
tt = {'IL', 'TR', 'PT'};

% % % ep_ind = ep_ind(1:18, :);

figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
for cond = 1:3
    nep = sum(ep_ind(:, cond));
% % %     carray_THVF = [logspace(0, 0.9, nep), logspace(0, 0.9, nep)]./10; % 0 is black line, 1 is white line, the darker the earlier line
% % %     carray = logspace(0, 0.9, nep)./10; % the darker the earlier line
    carray_THVF = [linspace(0, 0.9, nep), linspace(0, 0.9, nep)]; % the darker the earlier line
    carray = linspace(0, 0.9, nep); % the darker the earlier line
    
    subplot(2, 3, cond)
    p = polarplot(f_ang_th(ep_ind(:, cond), :)', f_mag_th(ep_ind(:, cond), :)', f_ang_vf(ep_ind(:, cond), :)', f_mag_vf(ep_ind(:, cond), :)', '--');
    for i = 1:length(p)
        p(i).Color = [carray_THVF(i), carray_THVF(i), carray_THVF(i)];
    end
    title(tt{cond})
    set(gca, 'FontSize', 18)
    
    subplot(2, 3, cond + 3)
    p = plot(tmp_lft_win_time(ep_ind(:, cond), :)', d_copy_thvf(ep_ind(:, cond), :)');
    for i = 1:length(p)
        p(i).Color = [carray(i), carray(i), carray(i)];
    end
    ylim([-50, 50])
% % %     vline(tch_time, ':k')
    xlim([-1, 2])
    vline(0, '--r', 'lift')
    ylabel('{\Delta}COPy_{TH-VF} (mm)');
end
mtit(sub_id)
%}

%% plot time profile I
%{
% % % ep_plot = 1:6;
% % % ep_plot = 7:12;
% % % ep_plot = 13:18;
% % % ep_plot = 19:24;
% % % ep_plot = 25:30;
% % % ep_plot = 31:36;
% % % ep_plot = 37:42;
% % % ep_plot = 43:48;
ep_plot = 49:54;
% % % ep_plot = 55:60;
% % % ep_plot = 61:66;
nep = length(ep_plot);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
id_plot = 1;
for ep = ep_plot
    subplot(3, nep, id_plot)
    plot(lft_win_time, f_mag_th(ep, :), lft_win_time, f_mag_vf(ep, :))
    ylim([0, 50])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('force (N)')
    end
    title(['t', num2str(data{ep, 1}), ' cond ', data{ep, 2}, ' ', num2str(data{ep, 4})])
    
    subplot(3, nep, nep + id_plot )
    plot(lft_win_time, rad2deg(f_ang_th(ep, :)), lft_win_time, rad2deg(f_ang_vf(ep, :)))
    ylim([-90, 270])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('angle ({\circ})')
        legend({'TH', 'VF'}, 'location', 'best')
    end
    
    subplot(3, nep, 2 * nep + id_plot)
    plot(lft_win_time, d_copy_thvf(ep, :));
    ylim([-30, 30])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('{\Delta}COPy_{TH-VF} (mm)');
    end
    xlabel('time (s)')
    id_plot = id_plot + 1;
end
suptitle(sub_id)
%}

%% plot time profile II
%{
ep_plot = 1:30; % 19:24; % 7:12; %1:6; % 7:12; % 13:18;
nep = length(ep_plot);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])

for ep = ep_plot
    subplot(5, 6, ep)
    plot(lft_win_time, f_mag_th(ep, :), lft_win_time, f_mag_vf(ep, :))
    ylim([0, 50])
    xlim([-1, 2])
    
    if ep == 1
        ylabel('force (N)')
        vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    else
        vline([0, tch_time(ep, 1)], {'--r', '--k'})
    end
    if ep <= 6
        title(['t', num2str(data{ep, 1}), ' cond ', data{ep, 2}, ' ', num2str(data{ep, 4})])
    end
end
suptitle(sub_id)
%}

%% time profile for Fn, dFy, dCOPy
%{
ep_plot = 1:6;
% % % ep_plot = 7:12;
% % % ep_plot = 13:18;
% % % ep_plot = 19:24;
% % % ep_plot = 25:30;
% % % ep_plot = 31:36;
% % % ep_plot = 37:42;
% % % ep_plot = 43:48;
% % % ep_plot = 49:54;
% % % ep_plot = 55:60;
% % % ep_plot = 61:66;
nep = length(ep_plot);
figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
id_plot = 1;
for ep = ep_plot
    subplot(3, nep, id_plot)
    plot(lft_win_time, f_n_th(ep, :), lft_win_time, -f_n_vf(ep, :))
    ylim([-1, 45])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('Fn (N)')
        legend({'TH', 'VF'}, 'location', 'best')
    end
    title(['t', num2str(data{ep, 1}), ' cond ', data{ep, 2}, ' ', num2str(data{ep, 4})])
    
    subplot(3, nep, nep + id_plot )
    plot(lft_win_time, d_fy_thvf(ep, :))
    ylim([-1, 15])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('{\Delta}Fy_{TH-VF} (N)');
    end
    
    subplot(3, nep, 2 * nep + id_plot)
    plot(lft_win_time, d_copy_thvf(ep, :));
    ylim([-20, 20])
    xlim([-1, 2])
    vline([0, tch_time(ep, 1)], {'--r', '--k'}, {'lift', 'touch'})
    if id_plot == 1
        ylabel('{\Delta}COPy_{TH-VF} (mm)');
    end
    xlabel('time (s)')
    id_plot = id_plot + 1;
end
suptitle(sub_id)
%}

%% for ploting dFy vs dCOPy and Fn vs dCOPy at lift onset
%{
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
[~, lft_onset] = min(abs(lft_win_time));

figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
ep_plot = 1:length(tmp_block); %5:length(tmp_block);
subplot(3, 3, 1)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), d_fy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), d_fy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), d_fy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), linespec)
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end          
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, d_fy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}
end
% % % legend('IL', 'TR', 'PT');
% sub-04: 18, 15, 12; sub-02: 13, 10, 7; sub-01: 18, 13, 8
text(-20, 18, 'IL', 'Color', 'r', 'FontSize', 18);
text(-20, 15, 'TR', 'Color', 'b', 'FontSize', 18);
text(-20, 12, 'PT', 'Color', 'k', 'FontSize', 18);

xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('{\Delta}Fy_{TH-VF} (N)');
vline(0, '--k')
hline(0, '--k')
hold off

subplot(3, 3, 2)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_n_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_n_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_n_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset), linespec)
    
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, f_n_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('Fn_{TH} (N)');
vline(0, '--k')
hold off

subplot(3, 3, 3)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_n_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_n_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('Fn_{VF} (N)');
vline(0, '--k')
hold off

subplot(3, 3, 5)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_mag_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_mag_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_mag_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset), linespec)
    
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, f_mag_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('F_{TH} (N)');
vline(0, '--k')
hold off

subplot(3, 3, 4)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), rad2deg(f_ang_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset)), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), rad2deg(f_ang_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset)), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), rad2deg(f_ang_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset)), linespec)
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, rad2deg(f_ang_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset)), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('angle_{TH} ({\circ})');
vline(0, '--k')
% % % hline(0, '--k')
hold off

subplot(3, 3, 8)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_mag_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_mag_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_mag_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), linespec)
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, f_mag_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}    
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('F_{VF} (N)');
vline(0, '--k')
hold off

subplot(3, 3, 7)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), rad2deg(f_ang_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset)), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), rad2deg(f_ang_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset)), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
%{    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), rad2deg(f_ang_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset)), linespec)
    if (ep < 4) || (ep > length(tmp_block) - 3)
        for i = 1:length(tmp_block{ep, 1})
            if i == 1
                tc = 'r';
            elseif i == length(tmp_block{ep, 1})
                tc = 'b';
            else
                tc = 'k';
            end
            text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset) + 0.3, rad2deg(f_ang_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset)), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tc, 'FontSize', 18)
        end
    end
%}
end
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('angle_{VF} ({\circ})');
vline(0, '--k')
% % % hline(0, '--k')
hold off

subplot(3, 3, 6)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_y_th([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_y_th([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
end
ylim([-10, 15])
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('Fy_{TH} (N)');
vline(0, '--k')
hold off

subplot(3, 3, 9)
hold on
for ep = ep_plot
    switch info_trial{tmp_block{ep, 1}(1, 1), 2}
        case 'IL'
            linespec = 'ro';
            tcc = 'r';
        case 'TR'
            linespec = 'bx';
            tcc = 'b';
        case 'PT'
            linespec = 'k^';
            tcc = 'k';
    end
    
    plot(d_copy_thvf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), f_y_vf([info_trial{tmp_block{ep, 1}, 1}], lft_onset), 'LineStyle', 'none', 'Marker', 'none')
    for i = 1:length(tmp_block{ep, 1})
        if (i == 1)
            fw = 'bold';
            fs = 12;
        elseif (i == length(tmp_block{ep, 1}))
            fw = 'bold';
            fs = 6;
        else
            fw = 'normal';
            fs = 3; %length(tmp_block{ep, 1}) - i;
        end

        text(d_copy_thvf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), f_y_vf([info_trial{tmp_block{ep, 1}(i), 1}], lft_onset), num2str([data{tmp_block{ep, 1}(i), 1}]), 'Color', tcc, 'FontWeight', fw, 'FontSize', fs)
    end
end
ylim([-10, 15])
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('Fy_{VF} (N)');
vline(0, '--k')
hold off

suptitle(All_info_trial{sub, 1})
%}