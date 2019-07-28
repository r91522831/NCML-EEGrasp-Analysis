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


%% calculate p-values
df = nsub - 1;
robust_tstat = mu ./ (sqrt(sig ./ nsub));

p = nan(nfreq, ntime, ncoeff, nelectrode);
for b = 1:ncoeff
    for e = 1:nelectrode
         p(:, :, b, e) = tcdf((robust_tstat(:, :, b, e)), df, 'upper'); %% change to correct dof ???? 
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

disp('here are some example plots')
%%

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