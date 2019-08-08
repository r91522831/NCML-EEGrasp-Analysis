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
coeff_name = {'b0', 'b1', 'b2', 'b3', 'b4', 'b5'};
ncoeff = length(coeff_name);

% Get data
dd = dir(fullfile(All_dirpath, '*LinearModel_coeff.mat'));
nsub = length(dd);
elec_dat = nan(nfreq, ntime, ncoeff, nsub, nelectrode);
for s = 1:nsub
    load(dd(s).name, 'Model_coeff_est');
    for e = 1:nelectrode
        elec_dat(:, :, :, s, e) = Model_coeff_est{e, 1};
    end
end
rmpath(All_dirpath)

%% TOO MANY FREQUENCIES - remove every other
% % % elec_dat = elec_dat(2:2:68, :, :, :, :);
% % % elec_dat = elec_dat(3:end, :, :, :, :);
nfreq = size(elec_dat, 1);

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
plot_preliminary('\alpha', 'C3', 'CP3', 9, 12, 450, 650, electrodes_name, freqz, timerstamps, mu, err_range); % alpha
legend('IL', 'TR', 'PT', 'Location', 'best')
subplot(2, 2, 3)
plot_preliminary('\beta', 'C1', 'C3', 20, 30, 450, 650, electrodes_name, freqz, timerstamps, mu, err_range); % beta
subplot(2, 2, 2)
plot_preliminary('\theta', 'FZ', 'FCz', 4, 8, 150, 350, electrodes_name, freqz, timerstamps, mu, err_range); % theta


%%
load(fullfile(fileparts(All_dirpath), 'preliminary_result.mat'));
load(fullfile(fileparts(All_dirpath), 'err_range.mat'));

% plot coefficient in electrode x freq animation
beta_0 = squeeze(mu(:, :, 1, :));

for e = 1:nelectrode
    for i = 1:nfreq
        for j = 1:ntime
            
        end
    end
end



%{
nf = 20;
nt = 50;
for e = 1%:nelectrode
    for i = 1:nf%nfreq
        for j = 1:nt%ntime
            % IL: beta_0 + beta_3 x Err; % TR: beta_1 + beta_4 x Err % PT: beta_2 + beta_5 x Err
            ersp_est_IL = mu(i, j, 1, e) + mu(i, j, 4, e) * err_range(:, j);
            ersp_est_TR = mu(i, j, 2, e) + mu(i, j, 5, e) * err_range(:, j);            
                
            subplot(nf, nt, (i - 1) * nt + j) %nfreq, ntime, i * ntime + j)
            plot(err_range(:, j) * 1000, ersp_est_IL, '-r', err_range(:, j) * 1000, ersp_est_TR, '-b');
            ylim([-1, 1])
            xlim([-5, 5])
            set(gca, 'xticklabels', [])
            set(gca, 'yticklabels', [])
%                 ylabel('estERSP change (%)')
%                 xlabel('err ({\circ})'
        end
    end
    
end
%}       




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