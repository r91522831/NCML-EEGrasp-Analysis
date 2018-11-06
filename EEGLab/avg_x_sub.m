close; clear; clc;

data_dir = uigetdir;
file_list = dir(fullfile(data_dir, '*.mat'));

sub_no = length(file_list);

data = nan(3, 6, sub_no);
data_abs_tf = nan(3, 6, sub_no);
for i = 1:sub_no
    load(fullfile(data_dir, file_list(i).name));
    
    data(:, :, i) = table2array(tf_power_complex);
    data_abs_tf(:, :, i) = table2array(tf_power);
end

%%
tf_x_sub_avg = abs(mean(data, 3));
tf_x_sub_std = std(data, 0, 3);

%% plot
figure(1)
set(0,'defaultAxesFontSize', 18)
subplot(2, 2, 1)
errorbar(tf_x_sub_avg(2, 2:end), tf_x_sub_std(2, 2:end)/sqrt(sub_no), 'o');
ylabel('Power change (%)')
xlim([0, 6])
xticks([1, 2, 3, 4, 5])
xticklabels({'IL', 'TR', 'PT_1', 'PT_2', 'PT_3'})
ylim([-2, 7])
title('\alpha', 'FontWeight', 'bold')
subplot(2, 2, 2)
errorbar(tf_x_sub_avg(1, 2:end), tf_x_sub_std(1, 2:end)/sqrt(sub_no), 'x');
ylabel('Power change (%)')
xlim([0, 6])
xticks([1, 2, 3, 4, 5])
xticklabels({'IL', 'TR', 'PT_1', 'PT_2', 'PT_3'})
ylim([-2, 7])
title('\theta', 'FontWeight', 'bold')
subplot(2, 2, 3)
errorbar(tf_x_sub_avg(3, 2:end), tf_x_sub_std(3, 2:end)/sqrt(sub_no), '*');
ylabel('Power change (%)')
xlim([0, 6])
xticks([1, 2, 3, 4, 5])
xticklabels({'IL', 'TR', 'PT_1', 'PT_2', 'PT_3'})
ylim([-2, 7])
title('\beta', 'FontWeight', 'bold')
