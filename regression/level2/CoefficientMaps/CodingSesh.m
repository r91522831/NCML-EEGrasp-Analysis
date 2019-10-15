close all; clear; clc;
All_dirpath = uigetdir();

addpath(All_dirpath);
% % % addpath(genpath(cd))
% % % cd('CoefficientMaps')
% Get constants
load('misc.mat')
% % % rmpath(genpath(cd))

timerstamps = tf_times{1}(:, :, 1);
freqz = squeeze(tf_freqs{1}(:, :, 1))';
% % % freqz = freqz(2:2:68);
% % % freqz = freqz(3:end);
ntime = length(timerstamps);
nfreq = length(freqz);
nelectrode = length(electrodes_name);
coeff_name = {'\beta_0', '\beta_1', '\beta_2', '\beta_3', '\beta_4', '\beta_5'};
ncoeff = length(coeff_name);

% Get data
dd = dir(fullfile(All_dirpath, '*LinearModel_coeff.mat'));

nsub = length(dd);
elec_dat = nan(nfreq, ntime, ncoeff, nsub, nelectrode);
ncond = 5; % IL1, IL2-19, TR1, TR2-19, and PT
roll_ang = cell(nsub, 2 + ncond);
for s = 1:nsub
    load(dd(s).name, 'Model_coeff_est', 'Model_Cond', 'Model_Roll');
    for e = 1:nelectrode
        elec_dat(:, :, :, s, e) = Model_coeff_est{e, 1};
    end
    % the Model_Roll angles are in radian
    roll_ang{s, 2} = [cellstr(Model_Cond), mat2cell(rad2deg(Model_Roll), ones(size(Model_Roll, 1), 1))];
    roll_ang{s, 1} = dd(s).name(1:6);
end
rmpath(All_dirpath)

% calculate roll angle for IL1, IL2-19, TR1, TR2-19, and PT
tmp = cell(1, ncond);
for s = 1:nsub
    if mod(str2double(roll_ang{s, 1}(end-1:end)), 2) ~= 1
        flip_side = 1;
    else
        flip_side = -1;
    end
    
    % IL1
    roll_ang{s, 2 + 1} = flip_side .* cell2mat(roll_ang{s, 2}(find(cell2mat(roll_ang{s, 2}(:, 1)) == '1', 1), 2));
    % IL2-19
    tmp_ind = find(cell2mat(roll_ang{s, 2}(:, 1)) == '1');
    roll_ang{s, 2 + 2} = flip_side .* cell2mat(roll_ang{s, 2}(tmp_ind(2:end, :), 2));
    % TR1
    roll_ang{s, 2 + 3} = flip_side .* cell2mat(roll_ang{s, 2}(find(cell2mat(roll_ang{s, 2}(:, 1)) == '2', 1), 2));
    % TR2-19
    tmp_ind = find(cell2mat(roll_ang{s, 2}(:, 1)) == '2');
    roll_ang{s, 2 + 4} = flip_side .* cell2mat(roll_ang{s, 2}(tmp_ind(2:end, :), 2));
    % PT
    roll_ang{s, 2 + 5} = flip_side .* cell2mat(roll_ang{s, 2}(cell2mat(roll_ang{s, 2}(:, 1)) == '3', 2));
end

sig_roll = nan(ncond, ntime);
mu_roll = nan(ncond, ntime);
for c = 1:ncond
    mu_roll(c, :) = mean(cell2mat(roll_ang(:, 2 + c)), 1);
    sig_roll(c, :) = var(cell2mat(roll_ang(:, 2 + c)), 0, 1);
end
save(fullfile(All_dirpath, 'roll.mat'), 'sig_roll', 'mu_roll', 'roll_ang')
% plot
%{
figure
h1 = shadedErrorBar(timerstamps, mu_roll(1, :), sqrt(sig_roll(1, :)), '-r', 1);
hold on
h2 = shadedErrorBar(timerstamps, mu_roll(2, :), sqrt(sig_roll(2, :)), '--m', 1);
h3 = shadedErrorBar(timerstamps, mu_roll(3, :), sqrt(sig_roll(3, :)), '-b', 1);
h4 = shadedErrorBar(timerstamps, mu_roll(4, :), sqrt(sig_roll(4, :)), '--c', 1);
h5 = shadedErrorBar(timerstamps, mu_roll(5, :), sqrt(sig_roll(5, :)), '--w', 1);
ylim([-20, 20])
xlim([-300, 1300])
vline(0, '--r')
xlabel('time (ms)')
ylabel('Roll angle ({\circ})')
hold off
lgnd = legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine, h5.mainLine], 'IL1', 'IL', 'TR1', 'TR', 'PT');
set(lgnd,'color','none');
set(gca, 'Color', [.8, .8, .8])
title('Object roll trajectories average across subjects')
%}

% delta: 1.5 ~ 4 Hz; theta: 4 ~ 8 Hz; alpha: 9 ~ 12 Hz; low beta: 13 ~ 19 Hz; high beta: 20 ~ 30 Hz; low gama: 30 ~ 35 Hz
rg_freq_band = { {'\delta', '1.5-4 Hz'}, {'\theta', '4-8 Hz'}, {'\alpha', '9-12 Hz'}, {'\beta_{low}', '13-19 Hz'}, {'\beta_{high}', '20-30 Hz'}, {'\gamma_{low}', '30-35 Hz'}; ...
                 find(freqz >=  2 & freqz <=  4), find(freqz >  4 & freqz <=  8), ...
                 find(freqz >  8   & freqz <= 13), find(freqz > 13 & freqz <= 20), ...
                 find(freqz > 20   & freqz <= 30), find(freqz > 30 & freqz <= 35) };
rg_time_win = { {'-50 to 150 ms'}, {'150 to 350 ms'}, {'350 to 450 ms'}, {'450 to 650 ms'}, {'650 to 850 ms'}; ...
                find(timerstamps >= -50 & timerstamps < 150), find(timerstamps >= 150 & timerstamps < 350), ...
                find(timerstamps >= 350 & timerstamps < 450), find(timerstamps >= 450 & timerstamps < 650), ...
                find(timerstamps >= 650 & timerstamps < 850) };
rg_elec_pick = { {'FZ and FCz'}, {'C3 and CP3'}, {'C1 and C3'}; ...
                find(strcmp('FZ', electrodes_name) | strcmp('FCz', electrodes_name)), ...
                find(strcmp('C3', electrodes_name) | strcmp('CP3', electrodes_name)), ...
                find(strcmp('C1', electrodes_name) | strcmp('C3', electrodes_name)) };

            
            
%% ========================================================================
% The two sections below are for estERSP plot.
% The estERSP is estimated at the individual subject level.
% The estERSP vs error is then average across subjects.
%% Calculate coefficients average across freq bands and time windows for each subject, electrode, and condition
nfreqband = length(rg_freq_band);
ntimewin = length(rg_time_win);
nelecpick = length(rg_elec_pick);
model_coeff_banded_wined = nan(nfreqband, ntimewin, ncoeff, nsub, nelecpick);
ticker = 1;
h = waitbar(0, 'calculating average coefficients.');
total_iter = nfreqband * ntimewin * ncoeff * nsub * nelecpick;
for s = 1:nsub
    for b = 1:ncoeff
        for e = 1:nelecpick
            for i = 1:nfreqband
                for j = 1:ntimewin
                    tmp = elec_dat(rg_freq_band{2, i}, rg_time_win{2, j}, b, s, rg_elec_pick{2, e});
                    model_coeff_banded_wined(i, j, b, s, e) = mean(tmp(:));
                    
                    ticker = ticker + 1;
                    progress_precent = 100 * ticker / total_iter;
                    waitbar(ticker / total_iter, h, sprintf('calculating average coefficients. %2.2f %%', progress_precent));
                end
            end
        end
    end
end
close(h)
save(fullfile(All_dirpath, 'result_freq_banded_time_wined.mat'), 'model_coeff_banded_wined', 'rg_freq_band', 'rg_time_win', 'coeff_name')

%% estimate ERSP for each subject with regard to max and min error
load(fullfile(fileparts(All_dirpath), 'err_range.mat'));
% err_range is in radian
ncond = 3;
estERSP = cell(nfreqband, ntimewin, ncond, nsub, nelecpick);
err = cell(ntimewin, 1);
for j = 1:ntimewin
    err{j, 1} = [max(err_range(1, rg_time_win{2, j})); min(err_range(2, rg_time_win{2, j}))];
    for s = 1:nsub
        for e = 1:nelecpick
            for i = 1:nfreqband
                for c = 1:ncond
                    switch c
                        case 1 % IL1 and IL2-19
                            estERSP{i, j, c, s, e} = elec_dat(i, j, 1, s, e) + elec_dat(i, j, 4, s, e) * err{j, 1};
                        case 2 % TR1 and TR2-19
                            estERSP{i, j, c, s, e} = elec_dat(i, j, 2, s, e) + elec_dat(i, j, 5, s, e) * err{j, 1};
                        case 3 % PT
                            estERSP{i, j, c, s, e} = elec_dat(i, j, 3, s, e) + elec_dat(i, j, 6, s, e) * err{j, 1};
                    end
                end
            end
        end
    end
end
% robust average ERSP across subject
sig_banded_wined = cell(nfreqband, ntimewin, ncond, nelecpick);
mu_banded_wined = cell(nfreqband, ntimewin, ncond, nelecpick);

for c = 1:ncond
    for e = 1:nelecpick
        for i = 1:nfreqband
            for j = 1:ntimewin
                tmp = [estERSP{i, j, c, :, e}];
                [tmp_sig_min, tmp_mu_min] = robustcov(tmp(1, :), 'Method', 'fmcd');
                [tmp_sig_max, tmp_mu_max] = robustcov(tmp(2, :), 'Method', 'fmcd');
                
                sig_banded_wined{i, j, c, e} = [tmp_sig_min; tmp_sig_max];
                mu_banded_wined{i, j, c, e} = [tmp_mu_min; tmp_mu_max];
            end
        end
    end
end
% plot
figure
subplot(2, 2, 1)
f = 3; t = 4; e = 2;
plot(rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 1, e}) - 1), '-ob', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 2, e}) - 1), '-sr', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 3, e}) - 1), '--dr');
ylabel('estERSP (%)')
xlabel({'err ({\circ})', [rg_freq_band{1, f}{1, 1}, ': ', rg_elec_pick{1, e}{1, 1}, ', freq: ', rg_freq_band{1, f}{1, 1}, ', time: ', rg_time_win{1, t}{1, 1}]})
legend('IL', 'TR', 'PT', 'Location', 'best')
subplot(2, 2, 3)
f = 5; t = 4; e = 3;
plot(rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 1, e}) - 1), '-ob', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 2, e}) - 1), '-sr', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 3, e}) - 1), '--dr');
ylabel('estERSP (%)')
xlabel({'err ({\circ})', [rg_freq_band{1, f}{1, 1}, ': ', rg_elec_pick{1, e}{1, 1}, ', freq: ', rg_freq_band{1, f}{1, 1}, ', time: ', rg_time_win{1, t}{1, 1}]})
subplot(2, 2, 2)
f = 2; t = 2; e = 1;
plot(rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 1, e}) - 1), '-ob', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 2, e}) - 1), '-sr', ...
     rad2deg(err{t, 1}), 100 * (db2pow(mu_banded_wined{f, t, 3, e}) - 1), '--dr');
ylabel('estERSP (%)')
xlabel({'err ({\circ})', [rg_freq_band{1, f}{1, 1}, ': ', rg_elec_pick{1, e}{1, 1}, ', freq: ', rg_freq_band{1, f}{1, 1}, ', time: ', rg_time_win{1, t}{1, 1}]})
%% ========================================================================




 





%% mild smoothing, need to account for in p-value
for s = 1:nsub
    for e = 1:nelectrode
        for b = 1:ncoeff
            elec_dat(:, :, b, s, e) = imgaussfilt(elec_dat(:, :, b, s, e), 1.28); % sigma = 1.28 => 2 voxel FWHM (full width at half maximum) smoothing
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FWHM = 2 * sqrt(2 * ln(2)) * sigma  %
% sigma = 0.84932 => FWHM = 2         %
% sigma = 1.274 => FWHM = 3           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate robust mean and sigma for freq bands
nfreqband = length(rg_freq_band);
sig_banded = nan(nfreqband, ntime, ncoeff, nelectrode);
mu_banded = nan(nfreqband, ntime, ncoeff, nelectrode);

ticker = 1;
h = waitbar(0, 'calculating robust means and covariances.');
total_iter = ncoeff * nelectrode * nfreqband * ntime;
for b = 1:ncoeff
    for e = 1:nelectrode
        for i = 1:nfreqband
            for j = 1:ntime
                tmp_band = squeeze(elec_dat(rg_freq_band{2, i}, j, b, :, e));
                [sig_banded(i, j, b, e), mu_banded(i, j, b, e)] = robustcov(tmp_band(:), 'Method', 'fmcd');
                ticker = ticker + 1;
                progress_precent = 100 * ticker / total_iter;
                waitbar(ticker / total_iter, h, sprintf('calculating robust means and covariances. %2.2f %%', progress_precent));
            end
        end
    end
end
close(h)
save(fullfile(All_dirpath, 'result_freq_banded.mat'), 'sig_banded', 'mu_banded', 'rg_freq_band', 'coeff_name')

%% Plot using topoplot
% /Users/yenhsunw/Dropbox (Personal)/Programming/Matlab/myLibrary/eeglab-develop/functions/sigprocfunc/topoplot.m
% % % All_dirpath = uigetdir();
load(fullfile(All_dirpath, 'misc.mat'));
load(fullfile(All_dirpath, 'result_freq_banded.mat'));
load(fullfile(All_dirpath, 'roll.mat'));
nfreqband = length(rg_freq_band);

coeff_trace = nan(ntime, ncoeff);
cmax = nan(1, ncoeff);
cmin = nan(1, ncoeff);
for b = 1:ncoeff
    coeff_trace(:, b) = mean(reshape(permute(squeeze(mu_banded(:, :, b, :)), [2, 1, 3]), 200, []), 2);
    tmp_mu = mu_banded(:, :, b, :);
    cmax(:, b) = max(tmp_mu(:));
    cmin(:, b) = min(tmp_mu(:));
end

% create the video writer with 1 fps
writerObj = VideoWriter(fullfile(All_dirpath, 'myVideo.avi'));
writerObj.FrameRate = 1;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
time_range = timerstamps >= -200 & timerstamps < 400;
figure(1)
for j = find(time_range)%1:ntime
    for b = 1:ncoeff
        for i = 1:nfreqband
            subplot(7, nfreqband, i + (b - 1) * nfreqband);
            % specify ('conv', 'on') to avoid extrapolation
            [tmp_h, ~, ~, ~, ~] = topoplot(mu_banded(i, j, b, :), chanlocs,'numcontour', 1, 'contourvals', z_score > critical, 'ccolor', 'w', 'maplimits', [cmin(1, b), cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
            if b == 1
                title(rg_freq_band{1, i})
            end
            if i == nfreqband
                colorbar;
            end
            if i == 1
                text(-1, 0, coeff_name{b});
            end
        end
    end
    subplot(7, nfreqband, [nfreqband * 6 + 1, nfreqband * 6 + 3])
    h1 = shadedErrorBar(timerstamps(time_range), mu_roll(1, time_range), sqrt(sig_roll(1, time_range)), '-r', 1);
    hold on
    h2 = shadedErrorBar(timerstamps(time_range), mu_roll(2, time_range), sqrt(sig_roll(2, time_range)), '-b', 1);
    h3 = shadedErrorBar(timerstamps(time_range), mu_roll(3, time_range), sqrt(sig_roll(3, time_range)), '-g', 1);
    h4 = shadedErrorBar(timerstamps(time_range), mu_roll(4, time_range), sqrt(sig_roll(4, time_range)), '-k', 1);
    hold off
    legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], 'IL', 'TR1', 'TR', 'PT', 'Location', 'eastoutside')
    vline(timerstamps(j), '-r')
    xlabel('time (ms)')
    ylabel('roll angle (\circ)')
    drawnow
    writeVideo(writerObj, getframe(gcf));
end
% close the writer object
close(writerObj);




%% ========================================================================
% The section below is not a good section since too many comparisons were performed!!!
%% plot one-sample t results and make a movie
%{
nfreqband = length(rg_freq_band);
% for one sample t-test
df = nsub - 1;
robust_tstat = mu_banded ./ (sqrt(sig_banded ./ nsub));

p_fbanded = nan(nfreqband, ntime, ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
         p_fbanded(:, :, b, e) = tcdf(robust_tstat(:, :, b, e), df, 'upper'); 
    end
end
% WHOLE FDR TxF within electrode
fdr = 0.15; % false discovery rate
%??? how to choose the desired false discovery rate? ??????????????????????
fp_all_fbanded = nan(nfreqband, ntime, ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
        p1_fbanded = p_fbanded(:, :, b, e);
        p1_fbanded = fdr_bh(p1_fbanded(:)', fdr);
        fp_all_fbanded(:, :, b, e) = reshape(p1_fbanded, nfreqband, size(p_fbanded, 2));
    end
end

% The significant mu array fdr-p-value > 0.15
fdr_mu_banded = fp_all_fbanded .* mu_banded;

cmax = nan(1, ncoeff);
cmin = nan(1, ncoeff);
for b = 1:ncoeff
    tmp_mu = fdr_mu_banded(:, :, b, :);
    cmax(:, b) = max(tmp_mu(:));
    cmin(:, b) = min(tmp_mu(:));
end

%% create the video writer with 1 fps
writerObj = VideoWriter(fullfile(All_dirpath, ['oneSample_pValue_fdr_point', num2str(fdr), '.avi']));
writerObj.FrameRate = 1;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
time_range = timerstamps >= -200 & timerstamps < 400;
figure('units','normalized','outerposition',[0 0 1 1])
for j = find(time_range)%1:ntime
    for b = 1:ncoeff
        for i = 1:nfreqband
            subplot(7, nfreqband, i + (b - 1) * nfreqband);
            % specify ('conv', 'on') to avoid extrapolation
            [tmp_h, ~, ~, ~, ~] = topoplot(fdr_mu_banded(i, j, b, :), chanlocs, 'maplimits', [cmin(1, b), cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
            if b == 1
                title(rg_freq_band{1, i})
            end
            if i == nfreqband
                colorbar;
            end
            if i == 1
                text(-1, 0, coeff_name{b});
            end
        end
    end
    subplot(7, nfreqband, [nfreqband * 6 + 1, nfreqband * 6 + 3])
    h1 = shadedErrorBar(timerstamps(time_range), mu_roll(1, time_range), sqrt(sig_roll(1, time_range)), '-r', 1);
    hold on
    h2 = shadedErrorBar(timerstamps(time_range), mu_roll(2, time_range), sqrt(sig_roll(2, time_range)), '-b', 1);
    h3 = shadedErrorBar(timerstamps(time_range), mu_roll(3, time_range), sqrt(sig_roll(3, time_range)), '-g', 1);
    h4 = shadedErrorBar(timerstamps(time_range), mu_roll(4, time_range), sqrt(sig_roll(4, time_range)), '-k', 1);
    hold off
    legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], 'IL', 'TR1', 'TR', 'PT', 'Location', 'eastoutside')
    vline(timerstamps(j), '-r')
    xlabel('time (ms)')
    ylabel('roll angle (\circ)')
    drawnow
    writeVideo(writerObj, getframe(gcf));
end
% close the writer object
close(writerObj);




%% for two sample t-test
nfreqband = length(rg_freq_band);
% H0: beta_3 = beta_4 (b = 4 and b = 5)
sig_banded_3 = squeeze(sig_banded(:, :, 4, :));
sig_banded_4 = squeeze(sig_banded(:, :, 5, :));
mu_banded_3 = squeeze(mu_banded(:, :, 4, :));
mu_banded_4 = squeeze(mu_banded(:, :, 5, :));
dmu_banded = mu_banded_3 - mu_banded_4;
robust_t2stat = dmu_banded ./ sqrt((sig_banded_3 .^2 ./ nsub) + (sig_banded_4 .^2 ./ nsub));
% The degrees of freedom nu associated with this variance estimate is approximated using the Welch–Satterthwaite equation
nu =   ((sig_banded_3 .^2 ./ nsub) + (sig_banded_4 .^2 ./ nsub)) .^2 ...
    ./ ((sig_banded_3 .^4 ./ (nsub ^2 * (nsub - 1))) + (sig_banded_4 .^4 ./ (nsub ^2 * (nsub - 1))));
p2_banded = nan(nfreqband, ntime, nelectrode);
for e = 1:nelectrode
     p2_banded(:, :, e) = tcdf(robust_t2stat(:, :, e), nu(:, :, e), 'upper'); 
end
% WHOLE FDR TxF within electrode
fdr = 0.35; % false discovery rate
%??? how to choose the desired false discovery rate? ??????????????????????
fp_all_fbanded_2 = nan(nfreqband, ntime, nelectrode);

for e = 1:nelectrode
    p1_fbanded = p2_banded(:, :, e);
    p1_fbanded = fdr_bh(p1_fbanded(:)', fdr);
    fp_all_fbanded_2(:, :, e) = reshape(p1_fbanded, nfreqband, size(p2_banded, 2));
end

% The significant mu array fdr-p-value > 0.15
fdr_mu_banded = fp_all_fbanded_2 .* mu_banded;



ncoeff = 1;
cmax = nan(1, ncoeff);
cmin = nan(1, ncoeff);
for b = 1:ncoeff
    tmp_mu = fdr_mu_banded(:, :, :);
    cmax(:, b) = max(tmp_mu(:));
    cmin(:, b) = min(tmp_mu(:));
end
%% create the video writer with 1 fps
writerObj = VideoWriter(fullfile(All_dirpath, 'twoSample_pValue.avi'));
writerObj.FrameRate = 1;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
time_range = timerstamps >= -200 & timerstamps < 400;
figure('units','normalized','outerposition',[0 0 1 1])
for j = find(time_range)%1:ntime
    for b = 1:ncoeff
        for i = 1:nfreqband
            subplot(7, nfreqband, i + (b - 1) * nfreqband);
            % specify ('conv', 'on') to avoid extrapolation
            [tmp_h, ~, ~, ~, ~] = topoplot(fdr_mu_banded(i, j, b, :), chanlocs, 'maplimits', [cmin(1, b), cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
            if b == 1
                title(rg_freq_band{1, i})
            end
            if i == nfreqband
                colorbar;
            end
            if i == 1
                text(-1, 0, coeff_name{b});
            end
        end
    end
    subplot(7, nfreqband, [nfreqband * 6 + 1, nfreqband * 6 + 3])
    h1 = shadedErrorBar(timerstamps(time_range), mu_roll(1, time_range), sqrt(sig_roll(1, time_range)), '-r', 1);
    hold on
    h2 = shadedErrorBar(timerstamps(time_range), mu_roll(2, time_range), sqrt(sig_roll(2, time_range)), '-b', 1);
    h3 = shadedErrorBar(timerstamps(time_range), mu_roll(3, time_range), sqrt(sig_roll(3, time_range)), '-g', 1);
    h4 = shadedErrorBar(timerstamps(time_range), mu_roll(4, time_range), sqrt(sig_roll(4, time_range)), '-k', 1);
    hold off
    legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], 'IL', 'TR1', 'TR', 'PT', 'Location', 'eastoutside')
    vline(timerstamps(j), '-r')
    xlabel('time (ms)')
    ylabel('roll angle (\circ)')
    drawnow
    writeVideo(writerObj, getframe(gcf));
end
% close the writer object
close(writerObj);

%}
%% ========================================================================












%% Calucluate robust mean and sigma
sig = nan(nfreq, ntime, ncoeff, nelectrode);
mu = nan(nfreq, ntime, ncoeff, nelectrode);
ticker = 1;
h = waitbar(0, 'calculating robust means and covariances.');
total_iter = ncoeff * nelectrode * nfreq * ntime;
for b = 1:ncoeff
    for e = 1:nelectrode
        for i = 1:nfreq
            for j = 1:ntime
                [sig(i, j, b, e), mu(i, j, b, e)] = robustcov(squeeze(elec_dat(i, j, b, :, e)), 'Method', 'fmcd');
                ticker = ticker + 1;
                progress_precent = 100 * ticker / total_iter;
                waitbar(ticker / total_iter, h, sprintf('calculating robust means and covariances. %2.2f %%', progress_precent));
            end
        end
    end
end
close(h)

save(fullfile(All_dirpath, 'preliminary_result.mat'), 'sig', 'mu')

%% Plot regression lines for each voxel (FxT) on ERSPxErr plan per electrode
load(fullfile(fileparts(All_dirpath), 'preliminary_result.mat'));
load(fullfile(fileparts(All_dirpath), 'err_range.mat'));

figure
subplot(2, 2, 1)
plot_preliminary('\alpha', 'C3', 'CP3', 9, 12, 450, 650, electrodes_name, freqz, timerstamps, mu, sig, err_range); % alpha
legend('IL', 'TR', 'PT', 'Location', 'best')
subplot(2, 2, 3)
plot_preliminary('\beta', 'C1', 'C3', 20, 30, 450, 650, electrodes_name, freqz, timerstamps, mu, sig, err_range); % beta
subplot(2, 2, 2)
plot_preliminary('\theta', 'FZ', 'FCz', 4, 8, 150, 350, electrodes_name, freqz, timerstamps, mu, sig, err_range); % theta




%% calculate p-values
%% for two sample t-test
% H0: beta_3 = beta_4 (b = 4 and b = 5)
sig_3 = squeeze(sig(:, :, 4, :));
sig_4 = squeeze(sig(:, :, 5, :));
mu_3 = squeeze(mu(:, :, 4, :));
mu_4 = squeeze(mu(:, :, 5, :));
dmu = mu_3 - mu_4;
robust_t2stat = dmu ./ sqrt((sig_3 .^2 ./ nsub) + (sig_4 .^2 ./ nsub));
% The degrees of freedom nu associated with this variance estimate is approximated using the Welch–Satterthwaite equation
nu =   ((sig_3 .^2 ./ nsub) + (sig_4 .^2 ./ nsub)) .^2 ...
    ./ ((sig_3 .^4 ./ (nsub ^2 * (nsub - 1))) + (sig_4 .^4 ./ (nsub ^2 * (nsub - 1))));
p2 = nan(nfreq, ntime, nelectrode);
for e = 1:nelectrode
     p2(:, :, e) = tcdf(robust_t2stat(:, :, e), nu(:, :, e), 'upper'); 
end

%% WHOLE FDR TxF within electrode
fp_all = nan(nfreq, ntime, nelectrode);
%??? how to choose the desired false discovery rate? ??????????????????????
fdr = 0.15; % false discovery rate
for e = 1:nelectrode
    p1 = p2(:, :, e);
    p1 = fdr_bh(p1(:)', fdr);
    fp_all(:, :, e) = reshape(p1, nfreq, size(p2, 2));
end

% % % b_name = {'IL', 'TR', 'PT', 'Err * IL', 'Err * TR', 'Err * PT'};

figure
for e = 1:nelectrode
    subplot(8, 8, e)
    contourf(timerstamps, freqz, fp_all(:, :, e) .* dmu(:, :, e), 'Linecolor', 'none');
    % % %         title(['cluster corrected mean ', coeff_name{b}, ' at electrode ', electrodes_name{e}]);
    title([electrodes_name{e}]);
end
suptitle('IL-TR regressor using whole electrode Time x frequency FDR')



















%% for one sample t-test
df = nsub - 1;
robust_tstat = mu ./ (sqrt(sig ./ nsub));

p = nan(nfreq, ntime, ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
         p(:, :, b, e) = tcdf(robust_tstat(:, :, b, e), df, 'upper'); 
    end
end

%% ???????? what does the 25 mean? ????????
%{
% ????????????????????????????????????????????????????????????????????????
nsub = 25;
df = nsub - 1;
robust_tstat = mu ./ (sqrt(sig ./ nsub));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtaining the p-value from the t stat and degree of freedom
% with Statistics Toolbox, use the p = tcdf(x, nu, 'upper') function
% otherwise use the lines below
% tdist2T = @(t,v) (1 - betainc(v / (v + t ^2), v / 2, 0.5));    % 2-tailed t-distribution
% tdist1T = @(t,v) 1 - (1 - tdist2T(t, v)) / 2;                  % 1-tailed t-distribution
% where ‘t’ is the t-statistic and ‘v’ are the degrees-of-freedom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:ncoeff
    for e = 1:nelectrode
% % %         p(:, :, b, e) = tdist2T(robust_tstat(:, :, i, j), df);
        p(:, :, b, e) = tcdf((robust_tstat(:, :, b, e)), df, 'upper'); %% change to correct dof
    end
end
% ????????????????????????????????????????????????????????????????????????
%}

%% Perform FDR across time in all electrodes per beta
time_idx = 35:104; % -100 ~ 600 ms? why select these time index? ?????????????????????
p = p(:, time_idx, :, :);
ntime = length(time_idx);

%% correct cluster sizes within electrode
cluster_pval = 0.05;
voxel_pval = 0.001;

p_indiv = nan(nfreq, ntime, ncoeff, nelectrode);
clustsizes = cell(ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
        clustinfo = bwconncomp( p(:, :, b, e) < voxel_pval );
        % identify clusters to remove
        clustsizes{b, e} = cellfun(@numel, clustinfo.PixelIdxList);
        clust_thresh = prctile( clustsizes{b, e}, 100 - cluster_pval * 100);
        whichclusters2remove = find(clustsizes{b, e} < clust_thresh );
        
        % remove clusters
        ptemp = p(:, :, b, e) < voxel_pval;
        for i = 1:length(whichclusters2remove)
            ptemp(clustinfo.PixelIdxList{whichclusters2remove(i)}) = 0;
        end
        p_indiv(:, :, b, e) = ptemp;
    end
end

%%
disp('here are some example plots')
for b = 1:ncoeff
    figure(b)
    for e = 1:nelectrode
        subplot(8, 8, e)
        contourf(timerstamps(time_idx), freqz, p_indiv(:, :, b, e) .* mu(:, time_idx, b, e), 'linecolor', 'none');
% % %         title(['cluster corrected mean ', coeff_name{b}, ' at electrode ', electrodes_name{e}]);
        title([electrodes_name{e}]);
    end
end


% % % figure;
% % % contourf(timerstamps(time_idx), freqz, p_indiv(:, :, 1, 2) .* mu(:, time_idx, 1, 2), 'linecolor', 'none');
% % % title('cluster corrected mean b0 at electrode Fcz');
% % % figure;
% % % contourf(timerstamps(time_idx), freqz, p_indiv(:, :, 4, 2) .* mu(:, time_idx, 4, 2), 'linecolor', 'none');
% % % title('cluster corrected mean Error*IL at electrode Fcz');



%{
%% I STOPPEd HERE
%% THIS IS A THRESHOLD OF CLUSTER SIZE USING CLUSTERS FROM ALL ELECTRODES
all_elec_thres = round(prctile([clustsizes{:}], 100 - cluster_pval * 100));

%%%%%% END USING CLUSTER CORRECTION APPROACH %%%%%%


%% You can just correct at particualr frequency bands using FDR
theta_pval = squeeze(mean(p(2:5, :, :, :), 1));
beta_pval = squeeze(mean(p(12:17, :, :, :), 1));
beta_high_pval = squeeze(mean(p(18:28, :, :, :), 1));
error_rate = 0.1; % The desired false discovery rate
for r = 1:ncoeff
    for e = 1:nelectrode
        p1 = theta_pval(:, r, e);
        p1 = fdr_bh(p1(:)', error_rate);
        theta_h1(:, r, e) = p1';
    end
end
figure;
imagesc(timerstamps(time_idx), 1:6, theta_h1(:, :, 2)');
ylabel('regressor')
xlabel('time (ms)')
title('theta (4-8 Hz) activity across regressors in electode Fcz')
disp('changing the 2 in theta_h1(:,:,2) changes the electrode')

%%
for r = 1:ncoeff
    for e = 1:nelectrode
        p1 = beta_pval(:, r, e);
            p1 = fdr_bh(p1(:)', error_rate);
            beta_h1(:, r, e) = p1';
    end
end
figure;
imagesc(timerstamps(time_idx), 1:6, beta_h1(:, :, 3)');
ylabel('regressor')
xlabel('time (ms)')
title('betalow(15-20 Hz) activity across regressors in electode C3')

%%
for r = 1:ncoeff
    for e = 1:nelectrode
        p1 = beta_high_pval(:, r, e);
            p1 = fdr_bh(p1(:)', error_rate);
            beta_high_h1(:, r, e) = p1';
    end
end
figure;
imagesc(timerstamps(time_idx), 1:6, beta_high_h1(:, :, 3)');
ylabel('regressor')
xlabel('time (ms)')
title('betahigh(20-30 Hz) activity across regressors in electode C3')
%}

%% WHOLE FDR TxF within electrode
fp_all = nan(nfreq, ntime, ncoeff, nelectrode);
for r = 1:ncoeff
    for e = 1:nelectrode
        p1 = p(:, :, r, e);
        p1 = fdr_bh(p1(:)', 0.15);
        fp_all(:, :, r, e) = reshape(p1, nfreq, size(p, 2));
    end
end

b_name = {'IL', 'TR', 'PT', 'Err * IL', 'Err * TR', 'Err * PT'};

for b = 1:ncoeff
    figure(b)
    for e = 1:nelectrode
        subplot(8, 8, e)
        contourf(timerstamps(time_idx), freqz, fp_all(:, :, b, e) .* mu(:, time_idx, b, e), 'Linecolor', 'none');
% % %         title(['cluster corrected mean ', coeff_name{b}, ' at electrode ', electrodes_name{e}]);
        title([electrodes_name{e}]);
    end
    suptitle([b_name{b}, ' regressor using whole electrode Time x frequency FDR'])
end
%% DIFFERENT APPROACH: set cluster threshold to k=20 and use FDR on all
%%% NOT COMPLETED

%% YOU COULD ALSO FIND THE PEAK WITHIN A TIME AND FREQUENCY WINDOW

%}