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

All_total_trial_n = 60;

nsub = length(All_selected_sub);
all_sub_mx_pRoll = nan(All_total_trial_n, 4, nsub);
cond_names = {'IL', 'TR', 'PT'};
nb_cond = length(cond_names);
mx_pRoll_cond = cell(nsub, nb_cond);

for sub = All_selected_sub%1:nsub
    clearvars -except sub filelist pathname all_sub_mx_pRoll nsub nb_cond mx_pRoll_cond cond_names All_*
    close all
    load(fullfile(pathname, filelist(sub).name));
    
    % add trial id to mx_onset and peak_roll
    tmp_filename = char({file_list(:).name});
    
    % get condition for trials
    tmp_cond_list = file_list;
    tmp_cond_id = cell(length(tmp_cond_list), 1);
    for j = 1:length(tmp_cond_list)
        tmp_ntrial = str2double(tmp_cond_list(j).name(7:9));
        if tmp_ntrial < 16
            tmp_cond_id{j, 1} = 'IL';
        elseif rem(tmp_ntrial - 15, 5) == 1
            tmp_cond_id{j, 1} = 'TR';
        else
            tmp_cond_id{j, 1} = 'PT';
        end
    end
    
    trial_id = str2double(cellstr(tmp_filename(:, 7:9)));
    tmp_mx_onset = mx_onset;
    tmp_peak_roll = peak_roll{:, 'peakRoll'}; % in degrees
    subID = str2double(filelist(sub).name(2:4));
    
    %{
    % aligned left and right
    if ~mod(subID, 2) % Odd number subject started IL with Left handle, while even number subject started with right handle!
        tmp_mx_onset = -mx_onset;
        tmp_peak_roll = -tmp_peak_roll;
    end
    %}
    
    tmp = nan(All_total_trial_n, 4);
    for i = 1:length(trial_id)
        tmp(trial_id(i), :) = [tmp_mx_onset(i), tmp_peak_roll(i), fy_onset(i), y_onset(i)];
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
ntrial = length(file_list);
session = cell(19, 2);
for i = 1:15
    session{1, 1} = 'I';
    session{1, 2}(i, :) = [i, avg_mx_pRoll(i, :), stde_mx_pRoll(i, :)];
end
% session(2, :) = {'T', [16, avg_mx_pRoll(16, :), stde_mx_pRoll(16, :)]};
j = 2;
tmp = [];
for i = 16:ntrial
    if rem(i, 5) == 2
        session(j, :) = {'T', tmp};
        j = j + 1;
        tmp = [];
    elseif i > 16 && rem(i, 5) == 1
        session(j, :) = {'P', tmp};
        j = j + 1;
        tmp = [];
    end
    
    tmp = [tmp; i, avg_mx_pRoll(i, :), stde_mx_pRoll(i, :)];
end
    
session(19, :) = {'P', tmp};

%% The expoential fit for the peak Roll
% % % All_beh_expfit = load(fullfile(pathname, ['behavior_pRoll_', num2str(nsub), 'subs_pRoll_fit.mat']));

%%
figure(1)
subplot 311
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
% % %     if strcmp(session{i, 1}, 'T')
% % %         errorbar(session{i, 2}(:, 1), -(session{i, 2}(:, 2)), session{i, 2}(:, 4), '-sb')
% % %     end
    errorbar(session{i, 2}(:, 1), session{i, 2}(:, 2), session{i, 2}(:, 4), line_spec)
end
hline(165, '--r')
hline(-165, '--r')
hold off
legend({'IL', 'TR', 'PT'}, 'Location', 'northwest')
% % % legend({'IL', 'TR', '-TR', 'PT'}, 'Location', 'northwest')
ylabel('peak Tcom (N-mm)')
% xlabel('trial')
% % % title(['error bars represent SE across ', num2str(length(All_selected_sub)), ' subjects'])
xticklabels([])
xlim([0, All_total_trial_n + 1])
set(gca, 'FontSize', 18)


subplot 312
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
% % %     if strcmp(session{i, 1}, 'T')
% % %         errorbar(session{i, 2}(:, 1), -(session{i, 2}(:, 2)), session{i, 2}(:, 4), '-sb')
% % %     end
    if strcmp(session{i, 1}, 'T')
        errorbar(session{i, 2}(:, 1), abs(session{i, 2}(:, 2) - 165), session{i, 2}(:, 4), line_spec)
    else
        errorbar(session{i, 2}(:, 1), abs(session{i, 2}(:, 2) + 165), session{i, 2}(:, 4), line_spec)
    end
end
hold off
legend({'IL', 'TR', 'PT'}, 'Location', 'northwest')
% % % legend({'IL', 'TR', '-TR', 'PT'}, 'Location', 'northwest')
ylabel('abs |peak Tcom_\epsilon| (N-mm)')
title('absolute Tcom error')
xticklabels([])
xlim([0, All_total_trial_n + 1])
set(gca, 'FontSize', 18)


subplot 313
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
% % %         yyaxis left
        h_err(i) = errorbar(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), session{i, 2}(:, 5), line_spec);
% % %         h_mirror(i) = errorbar(session{i, 2}(:, 1), session{i, 2}(:, 3), session{i, 2}(:, 5), '-oc');
    else
% % %         yyaxis left
        h(i) = shadedErrorBar(session{i, 2}(:, 1), abs(session{i, 2}(:, 3)), session{i, 2}(:, 5), line_spec);
        % exp fit for IL and PT
% % %         yyaxis right
% % %         pf(i) = plot(session{i, 2}(:, 1), All_beh_expfit.pRoll_dummy(session{i, 2}(:, 1), fit), line_spec_fit);
    end
end
% exp fit for TR
% % % yyaxis right
% % % pf(2) = plot(20:4:95, -All_beh_expfit.pRoll_dummy(20:4:end, 3), '-.b');

hold off
% % % yyaxis left
% % % ylim([-16, 16])
ylabel('peak roll ({\circ})')
% % % yyaxis right
% % % ylim([-1.4, 1.4])
% % % ylabel('exp fit dummy', 'rotation', 270, 'VerticalAlignment', 'bottom')
ylim([-10, 40])
% % % ylabel('peak roll ({\circ})')
legend([h(1).mainLine, h_err(16), h(17).mainLine], 'IL', 'abs(TR)', 'PT', 'Location', 'northwest')
% % % legend([h(1).mainLine, h_err(16), h_mirror(16), h(17).mainLine], 'IL', 'abs(TR)', 'TR', 'PT', 'Location', 'best')
% % % legend(pf(1:3), 'IL_{fit}', 'TR_{fit}', 'PT_{fit}', 'Location', 'southwest')
% % % hline(0, ':r')
%}
xlabel('trial')
xlim([0, All_total_trial_n + 1])


set(gca, 'FontSize', 18)
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])
savefig(fullfile(pathname, ['behavior_se_', num2str(nsub), 'subs']))