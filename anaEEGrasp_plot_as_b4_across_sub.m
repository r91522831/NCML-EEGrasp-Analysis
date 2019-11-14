%% Plot peak mx and peak roll around lift onset in all trials
close all; clearvars; clc

%% load aligned data
pathname = uigetdir;
filelist = dir(fullfile(pathname, '*_temp_result.mat'));

disp([num2cell((1:length(filelist))'), {filelist.name}']);
All_selected_sub = input('Which subject(s) to process? ');
if isempty(All_selected_sub)
    All_selected_sub = 1:length(filelist);
end

nsub = length(All_selected_sub);
all_sub_mx_pRoll = nan(95, 3, nsub);
cond_names = {'IL', 'TR', 'PT1', 'PT2', 'PT3'};
nb_cond = length(cond_names);
mx_pRoll_cond = cell(nsub, nb_cond);


for sub = All_selected_sub%1:nsub
    clearvars -except sub filelist pathname all_sub_mx_pRoll nsub nb_cond mx_pRoll_cond cond_names All_selected_sub
    close all
    load(fullfile(pathname, filelist(sub).name));
    
    % add trial id to mx_onset and peak_roll
    tmp_filename = char({file_list(:).name});
    
    % get condition for trials
    tmp_cond_list = file_list;
    tmp_cond_id = cell(length(tmp_cond_list), 1);
    for j = 1:length(tmp_cond_list)
        if strcmp(tmp_cond_list(j).name(11:12), 'PT')
            tmp_cond_id{j, 1} = tmp_cond_list(j).name([11:12, 18]);
        else
            tmp_cond_id{j, 1} = tmp_cond_list(j).name(11:12);
        end
    end
    
    trial_id = str2double(cellstr(tmp_filename(:, 7:9)));
    tmp_mx_onset = mx_onset;
    tmp_peak_roll = peak_roll{:, 'peakRoll'}; % in degrees
    subID = str2double(filelist(sub).name(2:4));
    % aligned left and right
    if ~mod(subID, 2) % Odd number subject started IL with Left handle, while even number subject started with right handle!
        tmp_mx_onset = -mx_onset;
        tmp_peak_roll = -tmp_peak_roll;
    end
    
    tmp = nan(95, 3);
    for i = 1:length(trial_id)
        tmp(trial_id(i), :) = [tmp_mx_onset(i), tmp_peak_roll(i), fy_onset(i)];
    end
    all_sub_mx_pRoll(:, :, sub) = tmp;
    
    for j = 1:nb_cond
        mx_pRoll_cond{sub, j} = tmp(strcmp(tmp_cond_id, cond_names{j}), :);
    end
end
avg_mx_pRoll = nanmean(all_sub_mx_pRoll(:, 1:2, :), 3);
stde_mx_pRoll = nanstd(all_sub_mx_pRoll(:, 1:2, :), 0, 3) ./ sqrt(nsub);

save(fullfile(pathname, ['behavior_pRoll_', num2str(nsub), 'subs']), 'all_sub_mx_pRoll', 'mx_pRoll_cond');

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

%% The expoential fit for the peak Roll
All_beh_expfit = load(fullfile(pathname, ['behavior_pRoll_', num2str(nsub), 'subs_pRoll_fit.mat']));

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
    if strcmp(session{i, 1}, 'T')
        errorbar(session{i, 2}(:, 1), -(session{i, 2}(:, 2)), session{i, 2}(:, 4), '-sc')
    end
    errorbar(session{i, 2}(:, 1), session{i, 2}(:, 2), session{i, 2}(:, 4), line_spec)
end
hold off
legend({'IL', '-TR', 'TR', 'PT'}, 'Location', 'northwest')
ylabel('peak Tcom (N-mm)')
% xlabel('trial')
title(['error bars represent SE across ', num2str(length(All_selected_sub)), ' subjects'])
xticklabels([])
xlim([0, 96])
set(gca, 'FontSize', 18)

subplot 212
hold on
for i = 1:length(session)
    switch session{i, 1}
        case 'I'
            line_spec = '-*r';
            line_spec_fit = '-.r';
            fit = 2;
        case 'T'
            line_spec = '-ob';
            line_spec_fit = '-.b';
            fit = 3;
        case 'P'
            line_spec = '-xk';
            line_spec_fit = '-.k';
            fit = 4;
        otherwise
    end
    if strcmp(session{i, 1}, 'T')
        yyaxis left
        h_err(i) = errorbar(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), session{i, 2}(:, 5), line_spec);
        h_mirror(i) = errorbar(session{i, 2}(:, 1), session{i, 2}(:, 3), session{i, 2}(:, 5), '-oc');
    else
        yyaxis left
        h(i) = shadedErrorBar(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), session{i, 2}(:, 5), line_spec);
        % exp fit for IL and PT
        yyaxis right
        pf(i) = plot(session{i, 2}(:, 1), All_beh_expfit.pRoll_dummy(session{i, 2}(:, 1), fit), line_spec_fit);
    end
end
% exp fit for TR
yyaxis right
pf(2) = plot(20:4:95, -All_beh_expfit.pRoll_dummy(20:4:end, 3), '-.b');

hold off
yyaxis left
ylim([-16, 16])
ylabel('peak roll ({\circ})')
yyaxis right
ylim([-1.4, 1.4])
ylabel('exp fit dummy', 'rotation', 270, 'VerticalAlignment', 'bottom')
% % % ylim([-20, 10])
% % % ylabel('peak roll ({\circ})')
% % % legend([h(1).mainLine, h_err(20), h_mirror(20), h(21).mainLine], 'IL', 'abs(TR)', 'TR', 'PT', 'Location', 'best')
legend(pf(1:3), 'IL_{fit}', 'TR_{fit}', 'PT_{fit}', 'Location', 'southwest')
xlabel('trial')
xlim([0, 96])
hline(0, ':r')
set(gca, 'FontSize', 18)
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])
savefig(fullfile(pathname, ['behavior_se_', num2str(nsub), 'subs']))