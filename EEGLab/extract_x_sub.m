clear; close all; clc;
for i = [4, 5, 6, 9]
    clearvars -except i
    load(['/Users/yen-hsunwu/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/forCPCA/S00', num2str(i), '.mat'])
    save(['S00', num2str(i), '_z_sub'], 'z_sub')
    clear z_sub
    save(['S00', num2str(i), '_h_sub'], 'h_sub')
    clear h_sub
    save(['S00', num2str(i), '_g_sub'], 'g_sub')
    clear g_sub
    save(['S00', num2str(i), '_tf_info'], 'tf_data', 'tf_ersp', 'tf_itc', 'tf_powbase', 'tf_times', 'tf_freqs')
end