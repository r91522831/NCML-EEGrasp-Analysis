close all; clear; clc;
All_dirpath = uigetdir();

addpath(All_dirpath);
% Get constants
load('misc.mat')

timerstamps = tf_times{1}(:, :, 1);
freqz = squeeze(tf_freqs{1}(:, :, 1))';
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the Model_Roll angles are in radian
    % rad2deg rad2deg rad2deg rad2deg rad2deg rad2deg rad2deg rad2deg
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
save(fullfile(All_dirpath, ['roll_', num2str(nsub),'.mat']), 'sig_roll', 'mu_roll', 'roll_ang')

%% plot
tmp_timewin = timerstamps >= -1000;
figure(ncoeff + 1);
h1 = shadedErrorBar(timerstamps(tmp_timewin) / 1000, mu_roll(1, tmp_timewin), sqrt(sig_roll(1, tmp_timewin)), '-r', 1);
hold on
h2 = shadedErrorBar(timerstamps(tmp_timewin) / 1000, mu_roll(2, tmp_timewin), sqrt(sig_roll(2, tmp_timewin)), '--m', 1);
h3 = shadedErrorBar(timerstamps(tmp_timewin) / 1000, mu_roll(3, tmp_timewin), sqrt(sig_roll(3, tmp_timewin)), '-b', 1);
h4 = shadedErrorBar(timerstamps(tmp_timewin) / 1000, mu_roll(4, tmp_timewin), sqrt(sig_roll(4, tmp_timewin)), '--c', 1);
h5 = shadedErrorBar(timerstamps(tmp_timewin) / 1000, mu_roll(5, tmp_timewin), sqrt(sig_roll(5, tmp_timewin)), '--y', 1);
% % % ylim([-20, 20])
% % % xlim([-300, 1300])
vline(0, '--r')
xlabel('time (s)')
ylabel('Roll angle ({\circ})')
hold off
lgnd = legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine, h5.mainLine], 'IL1', 'IL', 'TR1', 'TR', 'PT');
set(lgnd,'color','none');
set(gca, 'Color', [.8, .8, .8])
set(gca, 'FontSize', 24)
title(['Object roll trajectories average across ', num2str(nsub), ' subjects'])
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])
savefig(fullfile(All_dirpath, 'figures', ['roll_trajectory_', num2str(nsub)]))

%% delta: 2 ~ 4 Hz; theta: 4 ~ 8 Hz; alpha: 8 ~ 13 Hz; low beta: 13 ~ 20 Hz; high beta: 20 ~ 30 Hz; low gama: 30 ~ 35 Hz
rg_freq_band = { {'\delta', '2~4 Hz'}, {'\theta', '4~8 Hz'}, {'\alpha', '8~13 Hz'}, {'\beta_{low}', '13~20 Hz'}, {'\beta_{high}', '20~30 Hz'}, {'\gamma_{low}', '30~35 Hz'}; ...
                 find(freqz >=  2 & freqz <=  4), find(freqz >  4 & freqz <=  8), ...
                 find(freqz >   8 & freqz <= 13), find(freqz > 13 & freqz <= 20), ...
                 find(freqz >  20 & freqz <= 30), find(freqz > 30 & freqz <= 35) };
% % % rg_time_win = { {'-50 to 150 ms'}, {'150 to 350 ms'}, {'350 to 450 ms'}, {'450 to 650 ms'}, {'650 to 850 ms'}; ...
% % %                 find(timerstamps >= -50 & timerstamps < 150), find(timerstamps >= 150 & timerstamps < 350), ...
% % %                 find(timerstamps >= 350 & timerstamps < 450), find(timerstamps >= 450 & timerstamps < 650), ...
% % %                 find(timerstamps >= 650 & timerstamps < 850) };
rg_time_win = { {'-50~150', 'ms'}, {'150~250', 'ms'}, {'250~350', 'ms'}, {'350~450', 'ms'}, {'450~550', 'ms'}, {'550~650', 'ms'}; ...
                find(timerstamps >= -50 & timerstamps < 150), find(timerstamps >= 150 & timerstamps < 250), ...
                find(timerstamps >= 250 & timerstamps < 350), find(timerstamps >= 350 & timerstamps < 450), ...
                find(timerstamps >= 450 & timerstamps < 550), find(timerstamps >= 550 & timerstamps < 650) };
            
rg_elec_pick = { {'FZ and FCz'}, {'C3 and CP3'}, {'C1 and C3'}; ...
                find(strcmp('FZ', electrodes_name) | strcmp('FCz', electrodes_name)), ...
                find(strcmp('C3', electrodes_name) | strcmp('CP3', electrodes_name)), ...
                find(strcmp('C1', electrodes_name) | strcmp('C3', electrodes_name)) };
            
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
ntimeband = length(rg_time_win);
sig_banded = nan(nfreqband, ntimeband, ncoeff, nelectrode);
mu_banded = nan(nfreqband, ntimeband, ncoeff, nelectrode);

ticker = 1;
h = waitbar(0, 'calculating robust means and covariances.');
total_iter = ncoeff * nelectrode * nfreqband * ntimeband;
for b = 1:ncoeff
    for e = 1:nelectrode
        for i = 1:nfreqband
            for j = 1:ntimeband
                tmp_band = squeeze(elec_dat(rg_freq_band{2, i}, rg_time_win{2, j}, b, :, e));
                [sig_banded(i, j, b, e), mu_banded(i, j, b, e)] = robustcov(tmp_band(:), 'Method', 'fmcd');
                ticker = ticker + 1;
                progress_precent = 100 * ticker / total_iter;
                waitbar(ticker / total_iter, h, sprintf('calculating robust means and covariances. %2.2f %%', progress_precent));
            end
        end
    end
end
close(h)
save(fullfile(All_dirpath, ['result_ft_banded_', num2str(nsub), '.mat']), 'sig_banded', 'mu_banded', 'rg_freq_band', 'coeff_name')


%% for one sample t-test
df = nsub - 1;
robust_tstat = mu_banded ./ (sqrt(sig_banded ./ nsub));

p = nan(nfreqband, ntimeband, ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
         p(:, :, b, e) = tcdf(robust_tstat(:, :, b, e), df, 'upper'); 
    end
end


%% Plot using topoplot
% /Users/yenhsunw/Dropbox (Personal)/Programming/Matlab/myLibrary/eeglab-develop/functions/sigprocfunc/topoplot.m
% % % All_dirpath = uigetdir();
load(fullfile(All_dirpath, 'misc.mat'));
load(fullfile(All_dirpath, ['result_ft_banded_', num2str(nsub), '.mat']));
% % % load(fullfile(All_dirpath, 'result_freq_banded.mat'));
load(fullfile(All_dirpath, ['roll_', num2str(nsub), '.mat']));
nfreqband = length(rg_freq_band);
ntimeband = length(rg_time_win);

%
plotting_var = mu_banded;
% % % plotting_var = sqrt(sig_banded ./ nsub);
% % % plotting_var = robust_tstat;

% % % coeff_trace = nan(ntime, ncoeff);
cmax = nan(1, ncoeff);
cmin = nan(1, ncoeff);
for b = 1:ncoeff
% % %     coeff_trace(:, b) = mean(reshape(permute(squeeze(mu_banded(:, :, b, :)), [2, 1, 3]), 200, []), 2);
% % %     tmp_mu = mu_banded(:, :, b, :);
    tmp_mu = plotting_var(:, :, b, :);
    cmax(:, b) = max(tmp_mu(:));
    cmin(:, b) = min(tmp_mu(:));
end

if ~isfolder(fullfile(All_dirpath, 'figures'))
    mkdir(fullfile(All_dirpath, 'figures'));
end

critical = 0.05;
for b = 1:ncoeff
    figure(b)
    for j = 1:ntimeband
        for i = 1:nfreqband
            subplot(nfreqband, ntimeband + 1, j + (nfreqband - i) * (ntimeband + 1));
            % specify ('conv', 'on') to avoid extrapolation
            topoplot(plotting_var(i, j, b, :), chanlocs,'numcontour', 1, 'contourvals', p(i, j, b, :) < critical, 'ccolor', 'w', 'maplimits', [cmin(1, b), cmax(1, b)], 'electrodes', 'off', 'conv', 'on');
            if i == 1
                text(0, -0.8, rg_time_win{1, j}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
            end
            if j == 1
                text(-1.25, 0, rg_freq_band{1, i}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
            end
        end
    end
    subplot(nfreqband, ntimeband + 1, (ntimeband + 1) + ((1:nfreqband) - 1) .* (ntimeband + 1));
    surf([])
    tmp_c = max([abs(cmin(1, b)), abs(cmax(1, b))]);
	caxis([-tmp_c, tmp_c])
    set(gca, 'visible','off');
    colorbar;
    set(gca, 'FontSize', 16);
    set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])
    text(0, -0.1, ['coeff ', coeff_name{b}], 'FontWeight', 'bold', 'FontSize', 18);
    savefig(fullfile(All_dirpath, 'figures', ['coeff_beta_', num2str(b - 1), '_p_value_', num2str(nsub), 'subs']))
end

%% Plot avg across electrodes beta_0 + beta_3 * Roll for each conditions (IL, TR, PT) for Roll = -5 ~ 5
mu_beta = nan(nfreqband, ntimeband, ncoeff);
sig_beta = nan(nfreqband, ntimeband, ncoeff);
for i = 1:nfreqband
    for j = 1:ntimeband
        for b = 1:ncoeff
            tmp_beta = squeeze(mu_banded(i, j, b, :));
            [sig_beta(i, j, b), mu_beta(i, j, b) ]= robustcov(tmp_beta(:)); % beta_0
        end
    end
end

figure
for j = 1:ntimeband
    for i = 1:nfreqband
        subplot(nfreqband, ntimeband, j + (nfreqband - i) * ntimeband);
        x = -5:5;
        hold on
        % assume coeff_beta are not covary
% % %         h1 = plot(x, mu_beta(i, j, 1) + mu_beta(i, j, 4) * x, '-r');
% % %         h2 = plot(x, mu_beta(i, j, 2) + mu_beta(i, j, 5) * x, '-b');
% % %         h3 = plot(x, mu_beta(i, j, 3) + mu_beta(i, j, 6) * x, '-k');
        h1 = shadedErrorBar(x, mu_beta(i, j, 1) + mu_beta(i, j, 4) * x, sqrt(sig_beta(i, j, 1) + sig_beta(i, j, 4) * x.^2) ./ sqrt(nelectrode), '-r', 1);
        h2 = shadedErrorBar(x, mu_beta(i, j, 2) + mu_beta(i, j, 5) * x, sqrt(sig_beta(i, j, 2) + sig_beta(i, j, 5) * x.^2) ./ sqrt(nelectrode), '-b', 1);
        h3 = shadedErrorBar(x, mu_beta(i, j, 3) + mu_beta(i, j, 6) * x, sqrt(sig_beta(i, j, 3) + sig_beta(i, j, 6) * x.^2) ./ sqrt(nelectrode), '-k', 1);
        hold off
        if i == 1
            xlabel('roll ang ({\circ})');%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
        else
            xticklabels([])
            if i == nfreqband
                title(rg_time_win{1, j});%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
                if j == 1
% % %                     legend({'IL', 'TR', 'PT'});
                    legend([h1.mainLine, h2.mainLine, h3.mainLine], 'IL', 'TR', 'PT');
                end
            end
        end
        if j == 1
            ylabel(rg_freq_band{1, i});%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
        else
            yticklabels([])
        end
        ylim([-10, 10])
        set(gca, 'FontSize', 16);
    end
end
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);
% % % savefig(fullfile(All_dirpath, 'figures', 'est_power_across electrode'))
savefig(fullfile(All_dirpath, 'figures', ['est_power_across electrode_with se_', num2str(nsub), 'subs']))

%% plot regression for each electrode
% /Users/yenhsunw/Dropbox (Personal)/Programming/Matlab/myLibrary/eeglab-develop/functions/sigprocfunc/topoplot.m
% % % All_dirpath = uigetdir();
load(fullfile(All_dirpath, 'misc.mat'));
load(fullfile(All_dirpath, 'result_ft_banded.mat'));
% % % load(fullfile(All_dirpath, 'result_freq_banded.mat'));
load(fullfile(All_dirpath, 'roll.mat'));
nfreqband = length(rg_freq_band);
ntimeband = length(rg_time_win);

for e = 1:nelectrode
    figure
    for j = 1:ntimeband
        for i = 1:nfreqband
            subplot(nfreqband, ntimeband, j + (nfreqband - i) * ntimeband);
            x = -5:5;
            hold on
            % assume coeff_beta are not covary
            % % %         h1 = plot(x, mu_banded(i, j, 1, e) + mu_banded(i, j, 4, e) * x, '-r');
            % % %         h2 = plot(x, mu_banded(i, j, 2, e) + mu_banded(i, j, 5, e) * x, '-b');
            % % %         h3 = plot(x, mu_banded(i, j, 3, e) + mu_banded(i, j, 6, e) * x, '-k');
            h1 = shadedErrorBar(x, mu_banded(i, j, 1, e) + mu_banded(i, j, 4, e) * x, sqrt(sig_banded(i, j, 1, e) + sig_banded(i, j, 4, e) * x.^2) ./ sqrt(nsub), '-r', 1);
            h2 = shadedErrorBar(x, mu_banded(i, j, 2, e) + mu_banded(i, j, 5, e) * x, sqrt(sig_banded(i, j, 2, e) + sig_banded(i, j, 5, e) * x.^2) ./ sqrt(nsub), '-b', 1);
            h3 = shadedErrorBar(x, mu_banded(i, j, 3, e) + mu_banded(i, j, 6, e) * x, sqrt(sig_banded(i, j, 3, e) + sig_banded(i, j, 6, e) * x.^2) ./ sqrt(nsub), '-k', 1);
            hold off
            if i == 1
                xlabel('roll ang ({\circ})');%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
            else
                xticklabels([])
                if i == nfreqband
                    title(rg_time_win{1, j});%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
                    if j == 1
                        % % %                     legend({'IL', 'TR', 'PT'});
                        legend([h1.mainLine, h2.mainLine, h3.mainLine], 'IL', 'TR', 'PT');%, 'Location', 'northwestoutside');
                    end
                end
            end
            if j == 1
                ylabel(rg_freq_band{1, i});%, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
            else
                yticklabels([])
            end
            ylim([-10, 10])
            set(gca, 'FontSize', 16);
            vline(0, ':c')
            hline(0, ':c')
        end
    end
    suptitle(electrodes_name{e})
    set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);
    % % % savefig(fullfile(All_dirpath, 'figures', 'est_power_across electrode'))
    savefig(fullfile(All_dirpath, 'figures', ['est_power_', electrodes_name{e}, '_', num2str(nsub), 'subs']))
end






%{
% create the video writer with 1 fps
writerObj = VideoWriter(fullfile(All_dirpath, 'myVideo.avi'));
writerObj.FrameRate = 1;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
time_range = timerstamps >= -200 & timerstamps < 400;
critical = 0.05;
figure(1)
for j = find(time_range)%1:ntime
    for b = 1:ncoeff
        for i = 1:nfreqband
            subplot(7, nfreqband, i + (b - 1) * nfreqband);
            % specify ('conv', 'on') to avoid extrapolation
% % %             [tmp_h, ~, ~, ~, ~] = topoplot(mu_banded(i, j, b, :), chanlocs,'numcontour', 1, 'contourvals', z_score > critical, 'ccolor', 'w', 'maplimits', [cmin(1, b), cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
            [tmp_h, ~, ~, ~, ~] = topoplot(p(i, j, b, :), chanlocs,'numcontour', 1, 'contourvals', p(i, j, b, :) > critical, 'ccolor', 'w', 'maplimits', [cmin(1, b), cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
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













%{
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
%}
