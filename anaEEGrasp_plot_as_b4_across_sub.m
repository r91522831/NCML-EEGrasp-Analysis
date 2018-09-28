%% Plot peak mx and peak roll around lift onset in all trials
close all; clearvars; clc

%% load aligned data
pathname = uigetdir;
filelist = dir(fullfile(pathname, '*.mat'));

nsub = length(filelist);
all_sub_mx_pRoll = nan(95, 2, nsub);
for sub = 1:nsub
    clearvars -except sub filelist pathname all_sub_mx_pRoll nsub
    close all
    load(fullfile(pathname, filelist(sub).name));
    
    % add trial id to mx_onset and peak_roll
    tmp_filename = char({file_list(:).name});
    trial_id = str2num(tmp_filename(:, 7:9));
    tmp_mx_onset = mx_onset;
    tmp_peak_roll = peak_roll{:, 'peakRoll'};
    if ~mod(sub, 2)
        tmp_mx_onset = -mx_onset;
        tmp_peak_roll = -tmp_peak_roll;
    end
    if length(trial_id) ~= 95
        tmp = nan(95, 2);
        tmp_forward = 1;
        for i = 1:95
            if tmp_forward > length(trial_id)
                break;
            end
            if i ~= trial_id(tmp_forward)
                continue;
            end
            tmp(i, :) = [tmp_mx_onset(tmp_forward), tmp_peak_roll(tmp_forward)];
            tmp_forward = tmp_forward + 1;
        end
    else
        tmp = [tmp_mx_onset, tmp_peak_roll];
    end
    all_sub_mx_pRoll(:, :, sub) = tmp;
end
avg_mx_pRoll = nanmean(all_sub_mx_pRoll, 3);
stde_mx_pRoll = nanstd(all_sub_mx_pRoll, 0, 3);% ./ sqrt(nsub);
%%
tmp = [str2double(file_list(1).name(7:9)), avg_mx_pRoll(1, :), stde_mx_pRoll(1, :)];
ntrial = length(file_list);
session = cell(1, 2);
j = 1;
for i = 2:ntrial
    if (file_list(i).name(:, 11) ~= file_list(i - 1).name(:, 11))
        session(j, :) = {file_list(i - 1).name(:, 11), tmp};
        j = j + 1;
        tmp = [];
    end
    
    tmp = [tmp; str2double(file_list(i).name(7:9)), avg_mx_pRoll(i, :), stde_mx_pRoll(i, :)];
end
session(j, :) = {file_list(end).name(:, 11), tmp};

%%
figure(1)
subplot 211
hold on
for i = 1:length(session)
    switch session{i, 1}
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    errorbar(session{i, 2}(:, 1), session{i, 2}(:, 2), session{i, 2}(:, 4), line_spec)
end
hold off
legend({'IL', 'TR', 'PT'}, 'Location', 'southwest')
ylabel('Tcom (N-mm)')
xlabel('trial')
xlim([0, 96])
subplot 212
hold on
for i = 1:length(session)
    switch session{i, 1}
        case 'I'
            line_spec = '-*r';
        case 'T'
            line_spec = '-ob';
        case 'P'
            line_spec = '-xk';
        otherwise
    end
    if strcmp(session{i, 1}, 'T')
        errorbar(session{i, 2}(:, 1), session{i, 2}(:, 3), session{i, 2}(:, 5), line_spec)
    else
        shadedErrorBar(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), abs(session{i, 2}(:, 5)), line_spec)
    end
end
hold off
ylim([-10, 21])
ylabel('absolute peak roll ({\circ})')
xlabel('trial')
xlim([0, 96])
mtit('error bars are SD')
savefig(fullfile(pathname, 'behavior_sd'))