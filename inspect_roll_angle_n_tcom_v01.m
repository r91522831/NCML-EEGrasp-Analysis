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
    info_trial = cell(nep, 9);
    for ep = 1:length(beh.file_list)
        info_trial{ep, 1} = str2double(beh.file_list(ep).name(7:9)); % trial id
        info_trial{ep, 2} = beh.file_list(ep).name(11:12); % condition: IL, TR, PT, TR
        info_trial{ep, 3} = str2double(beh.file_list(ep).name(23:25)); % added weight in gram
        info_trial{ep, 4} = str2double(beh.file_list(ep).name(16:18)); % trial number in block
        
        deltaCOPy = beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'COPy'} - beh.finger_V{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'COPy'};
        info_trial{ep, 7} = deltaCOPy;
        info_trial{ep, 8} = sqrt(beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'}.^2 + beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fz'}.^2);
        info_trial{ep, 9} = beh.finger_Th{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'} - beh.finger_V{ep, 1}{beh.info_onset_time{ep, 'lft_ind'}, 'fy'};
    end
    
    info_trial(:, 5) = mat2cell([beh.peak_roll{:, 'peakRoll'}], ones(nep, 1));
    info_trial(:, 6) = mat2cell([beh.peak_mx{:, 'peakMx'}], ones(nep, 1));
    
    All_info_trial{All_sub_i, 1} = ['sub-', sub_id];
    All_info_trial{All_sub_i, 2} = info_trial;
    All_info_trial{All_sub_i, 3} = diff(beh.info_time_trigger{1, 1}(1:2, 1)); % dt in ms
    All_sub_i = All_sub_i + 1;
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
tmp_block{b, 1} = tmp_b:n_trial;

%% plot individual trials
f = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1], 'DefaultAxesFontSize', 18);
subplot 411
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

%{
% for 1st pilot
ylim([-700, 700])
session = [32.5, 64.5, 96.5, 128.5];
text([16, 48, 80, 112, 144], [700, 700, 700, 700, 700], {'0g', '300g', '0g', '150g', '0g'}, 'FontSize', 18);
plot([32.5, 64.5, 64.5, 96.5, 96.5, 128.5, 128.5, 160], [564, 564, 395, 395, 497, 497, 395, 395], 'r--')
plot([32.5, 64.5, 64.5, 96.5, 96.5, 128.5, 128.5, 160], -1 * [564, 564, 395, 395, 497, 497, 395, 395], 'r--')
vline(session, 'k--');
%}

hline([-395, 395], 'r--')
hold off
ylabel('Tcom (Nmm)');
title(All_info_trial{sub, 1});

subplot 412
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
    plot([info_trial{tmp_block{ep, 1}, 1}], [info_trial{tmp_block{ep, 1}, 7}], linespec);
end
xlim([0, n_trial + 1])
ylim([-40, 40])
hline(0)

% for 1st pilot
% % % vline(session, 'k--');

hold off
ylabel('{\Delta}COPy_{TH-VF} (mm)');

subplot 413
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
    plot([info_trial{tmp_block{ep, 1}, 1}], [info_trial{tmp_block{ep, 1}, 9}], linespec);
end
xlim([0, n_trial + 1])
hline(0)

% for 1st pilot
% % % vline(session, 'k--');

hold off
ylabel('{\Delta}Fy_{TH-VF} (N)');

subplot 414
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

% for 1st pilot
% % % vline(session, 'k--');

ylabel('peak roll ({\circ})');
xlabel('trial ID');


%% average across all blocks
[nrep, ~, id_rep]= unique(cellfun(@(x) length(x), tmp_block));

avg_ILPT = cell(length(nrep), 1);
avg_TR = cell(length(nrep), 1);
std_ILPT = avg_ILPT; std_TR = avg_TR;
for r = 1:length(nrep)
    data_ILPT = [];
    data_TR = [];
    ind_b = find((id_rep == r));
    for b = 1:length(ind_b)
        if mod(b, 2)
            data_ILPT = cat(3, data_ILPT, info_trial(tmp_block{ind_b(b), 1}, [5:7, 9]));
        else
            data_TR = cat(3, data_TR, info_trial(tmp_block{ind_b(b), 1}, [5:7, 9]));
        end
    end
    avg_ILPT{r} = mean(cell2mat(data_ILPT), 3);
    std_ILPT{r} = std(cell2mat(data_ILPT), [], 3);
    avg_TR{r} = mean(cell2mat(data_TR), 3);
    std_TR{r} = std(cell2mat(data_TR), [], 3);
end

%% plot avg across all blocks
f1 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1], 'DefaultAxesFontSize', 18);
subplot 411
hold on;
id_t = 1;
for r = length(nrep):-1:1
    nr = nrep(r);
    id_trial = id_t:(id_t + nr - 1);
% % %     plot(id_trial, avg_ILPT{r, 1}(:, 2), '-ro', id_trial, avg_TR{r, 1}(:, 2), '-bx');
    errorbar(id_trial, avg_ILPT{r, 1}(:, 2), std_ILPT{r, 1}(:, 2), '-ro');
    errorbar(id_trial, avg_TR{r, 1}(:, 2), std_TR{r, 1}(:, 2), '-bx');
    id_t = id_t + nr + 1;
end
ylim([-600, 600])
xlim([0, id_t - 1]);
hline([-395, 395], 'r--')
hold off
ylabel('Tcom (Nmm)');
legend({'R', 'L'}, 'location', 'best')
title(All_info_trial{sub, 1});

subplot 412
hold on;
id_t = 1;
for r = length(nrep):-1:1
    nr = nrep(r);
    id_trial = id_t:(id_t + nr - 1);
% % %     plot(id_trial, avg_ILPT{r, 1}(:, 3), '-ro', id_trial, avg_TR{r, 1}(:, 3), '-bx');
    errorbar(id_trial, avg_ILPT{r, 1}(:, 3), std_ILPT{r, 1}(:, 3), '-ro');
    errorbar(id_trial, avg_TR{r, 1}(:, 3), std_TR{r, 1}(:, 3), '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
ylim([-40, 40])
hline(0)
hold off
ylabel('{\Delta}COPy_{TH-VF} (mm)');

subplot 413
hold on;
id_t = 1;
for r = length(nrep):-1:1
    nr = nrep(r);
    id_trial = id_t:(id_t + nr - 1);
% % %     plot(id_trial, avg_ILPT{r, 1}(:, 4), '-ro', id_trial, avg_TR{r, 1}(:, 4), '-bx');
    errorbar(id_trial, avg_ILPT{r, 1}(:, 4), std_ILPT{r, 1}(:, 4), '-ro');
    errorbar(id_trial, avg_TR{r, 1}(:, 4), std_TR{r, 1}(:, 4), '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
hline(0)
hold off
ylabel('{\Delta}Fy_{TH-VF} (N)');

subplot 414
hold on;
id_t = 1;
for r = length(nrep):-1:1
    nr = nrep(r);
    id_trial = id_t:(id_t + nr - 1);
% % %     plot(id_trial, avg_ILPT{r, 1}(:, 1), '-ro', id_trial, avg_TR{r, 1}(:, 1), '-bx');
    errorbar(id_trial, abs(avg_ILPT{r, 1}(:, 1)), std_ILPT{r, 1}(:, 1), '-ro');
    errorbar(id_trial, abs(avg_TR{r, 1}(:, 1)), std_TR{r, 1}(:, 1), '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
ylim([0, 18])
ylabel('peak roll ({\circ})');
xlabel('trial ID');

%% average across several blocks
tmp_group = {1, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 4, 5, 6, 5, 6, 5, 6, 5, 6, 5}';
tmp_block(1:21, 2) = tmp_group;
tmp_block(22:42, 2) = tmp_group;
tmp_block(43:63, 2) = tmp_group;

data_21 = cell(6, 1);
for e = 1:21
    tmp_d = data_21{tmp_block{e, 2}, 1};
    data_21{tmp_block{e, 2}, 1} = cat(3, tmp_d, info_trial(tmp_block{e, 1}, [5:7, 9]));
end
data_42 = cell(6, 1);
for e = 22:42
    tmp_d = data_42{tmp_block{e, 2}, 1};
    data_42{tmp_block{e, 2}, 1} = cat(3, tmp_d, info_trial(tmp_block{e, 1}, [5:7, 9]));
end
data_63 = cell(6, 1);
for e = 43:63
    tmp_d = data_63{tmp_block{e, 2}, 1};
    data_63{tmp_block{e, 2}, 1} = cat(3, tmp_d, info_trial(tmp_block{e, 1}, [5:7, 9]));
end
data_all = [data_21; data_42; data_63];

avg_all = cell(length(data_all), 1);
std_all = cell(length(data_all), 1);
for i = 1:length(data_all)
    avg_all{i, 1} = mean(cell2mat(data_all{i, 1}), 3);
    std_all{i, 1} = std(cell2mat(data_all{i, 1}), [], 3);
end

%% plot avg across all blocks
f2 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1], 'DefaultAxesFontSize', 18);
subplot 411
hold on;
id_t = 1;
for i = 1:2:length(avg_all)
    avg_ILPT = avg_all{i, 1}(:, 2);
    std_ILPT = std_all{i, 1}(:, 2);
    avg_TR = avg_all{i + 1, 1}(:, 2);
    std_TR = std_all{i + 1, 1}(:, 2);
    nr = length(avg_ILPT);
    id_trial = id_t:(id_t + nr - 1);
    
% % %     plot(id_trial, avg_ILPT, '-ro', id_trial, avg_TR, '-bx');
    errorbar(id_trial, avg_ILPT, std_ILPT, '-ro');
    errorbar(id_trial, avg_TR, std_TR, '-bx');
    id_t = id_t + nr + 1;
end
ylim([-600, 600])
xlim([0, id_t - 1]);
hline([-395, 395], 'r--')
hold off
ylabel('Tcom (Nmm)');
legend({'R', 'L'}, 'location', 'best')
title(All_info_trial{sub, 1});

subplot 412
hold on;
id_t = 1;
for i = 1:2:length(avg_all)
    avg_ILPT = avg_all{i, 1}(:, 3);
    std_ILPT = std_all{i, 1}(:, 3);
    avg_TR = avg_all{i + 1, 1}(:, 3);
    std_TR = std_all{i + 1, 1}(:, 3);
    nr = length(avg_ILPT);
    id_trial = id_t:(id_t + nr - 1);
    
% % %     plot(id_trial, avg_ILPT, '-ro', id_trial, avg_TR, '-bx');
    errorbar(id_trial, avg_ILPT, std_ILPT, '-ro');
    errorbar(id_trial, avg_TR, std_TR, '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
ylim([-40, 40])
hline(0)
hold off
ylabel('{\Delta}COPy_{TH-VF} (mm)');

subplot 413
hold on;
id_t = 1;
for i = 1:2:length(avg_all)
    avg_ILPT = avg_all{i, 1}(:, 4);
    std_ILPT = std_all{i, 1}(:, 4);
    avg_TR = avg_all{i + 1, 1}(:, 4);
    std_TR = std_all{i + 1, 1}(:, 4);
    nr = length(avg_ILPT);
    id_trial = id_t:(id_t + nr - 1);
    
% % %     plot(id_trial, avg_ILPT, '-ro', id_trial, avg_TR, '-bx');
    errorbar(id_trial, avg_ILPT, std_ILPT, '-ro');
    errorbar(id_trial, avg_TR, std_TR, '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
hline(0)
hold off
ylabel('{\Delta}Fy_{TH-VF} (N)');

subplot 414
hold on;
id_t = 1;
for i = 1:2:length(avg_all)
    avg_ILPT = avg_all{i, 1}(:, 1);
    std_ILPT = std_all{i, 1}(:, 1);
    avg_TR = avg_all{i + 1, 1}(:, 1);
    std_TR = std_all{i + 1, 1}(:, 1);
    nr = length(avg_ILPT);
    id_trial = id_t:(id_t + nr - 1);
    
% % %     plot(id_trial, abs(avg_ILPT), '-ro', id_trial, abs(avg_TR), '-bx');
    errorbar(id_trial, abs(avg_ILPT), std_ILPT, '-ro');
    errorbar(id_trial, abs(avg_TR), std_TR, '-bx');
    id_t = id_t + nr + 1;
end
xlim([0, id_t - 1]);
ylim([0, 18])
ylabel('peak roll ({\circ})');
xlabel('trial ID');
