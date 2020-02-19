%%
sub_id = EEG.filename(1:6);
% electrode_id = find(strcmpi({EEG.chanlocs.labels}, 'CP4'));
cond = {'IL', 'TR', 'PT1'};
smooth = 2;

for electrode_i = 1:length({EEG.chanlocs.labels})
    % get color limit
    data = cell(length(cond), 1);
    for i = 1:length(cond)
        subplot(3, 1, i)
        data{i} = pop_erpimage(EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}], smooth, 1, {[cond{i}, '_onset']}, [], 'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', 'NoShow', 'on');
        get(gca, 'CLim')
    end
    data = cell2mat(data);
    clim = max(max(max(data), abs(min(min(data)))));
    colorlimit = [-clim, clim];
    % plot
    h = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    for i = 1:length(cond)
        subplot(3, 1, i)
        pop_erpimage(EEG, 1, [electrode_i], [[]], [EEG.chanlocs(electrode_i).labels, ' ', cond{i}], smooth, 1, {[cond{i}, '_onset']}, [], 'condID', 'yerplabel', '\muV', 'limits', [-2000 3000 NaN NaN NaN NaN NaN NaN], 'cbar', 'on', 'caxis', colorlimit);
    end
    
    electrode_label = EEG.chanlocs(electrode_i).labels;
    figdir = '/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/for Feb18 2020';
    figname = ['erpimage_', electrode_label, '_', sub_id];
    savefig(h, fullfile(figdir, figname))
    saveas(h, fullfile(figdir, [figname, '.png']))
    close(h)
end


%% for different frequency bands
% % %     tf = load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/eeglab/003 StartOver/linear/RAW/sub-XX_timefreq.mat');
tf = load(fullfile(All_path, All_filelist(All_i).name));

tf_ersp = tf.tf_ersp;
timerstamps = tf.tf_times';
freqz = tf.tf_freqs';
electrodes = tf.tf_ersp.Properties.RowNames;
nb_epoch = size(tf.tf_ersp{1, 1}{:}, 3);

ntime = length(timerstamps);
nfreq = length(freqz);
nelectrode = length(electrodes);

% % % rg_time = find(timerstamps >= 400 & timerstamps < 600); % 400 to 600 ms after lift onset

% theta band: 4 to 8 Hz; low beta: 13 ~ 19 Hz; high beta: 20 ~ 30 Hz
rg_freq = {find(freqz > 4 & freqz <= 8), find(freqz > 13 & freqz <= 20), find(freqz > 20 & freqz <= 30) };
rg_freq_name = {'\theta', 'low\beta', 'high\beta'};
nfreqband = length(rg_freq);

power_mean_fb = cell(1, nfreqband);
for i = 1:nelectrode
    for fb = 1:nfreqband
        power_mean_fb{1, fb}(i, :, :) = mean(tf_ersp{i, 1}{:}(rg_freq{fb}, :, :), 1);
    end
    
% % %     for t = 1:ntime
% % %         for ep = 1:nb_epoch
% % %             [~, power_theta{i, 1}(t, ep)] = robustcov(tf_ersp{i, 1}{:}(rg_freq{1}, t, ep));
% % %         end
% % %     end
end

% plot theta
originEEG = EEG;
EEG.data = power_mean_fb{:, 1};
EEG.times = timerstamps;




