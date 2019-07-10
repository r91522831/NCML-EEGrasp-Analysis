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
freqz = freqz(2:2:68);
freqz = freqz(3:end);

% Get data
dd = dir(fullfile(All_dirpath, '*LinearModel_coeff.mat'));
nfreq = 68; ntime = 200; ncoeff = 6; nsub = length(dd); nelectrode = 5;
elec_dat = nan(nfreq, ntime, ncoeff, nsub, nelectrode);
for i = 1:nsub
    load(dd(i).name, 'Model_coeff_est');
    for e = 1:nelectrode
        elec_dat(:, :, :, i, e) = Model_coeff_est{e, 1};
    end
end
rmpath(All_dirpath)

%% TOO MANY FREQUENCIES - remove every other
elec_dat = elec_dat(2:2:68, :, :, :, :);
elec_dat = elec_dat(3:end, :, :, :, :);

%% mild smoothing, need to account for in p-value
for i = 1:length(dd)
    for e = 1:nelectrode
        for s = 1:ncoeff
            elec_dat(:, :, s, i, e) = imgaussfilt(elec_dat(:, :, s, i, e), 1.28); % ? = 1.28 => 2 voxel FWHM (full width at half maximum) smoothing
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FWHM = 2 * sqrt(2 * ln(2)) * ?  %
% ? = 0.84932 => FWHM = 2         %
% ? = 1.274 => FWHM = 3           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calucluate robust mean and sigma
ticker = 1;
h = waitbar(0, 'calculating robust means and covariances.');
for b = 1:ncoeff
    for r = 1:nelectrode
        for i = 1:nfreq
            for j = 1:ntime
                [sig(i, j, b, r), mu(i, j, b, r), ~, ~] = robustcov(squeeze(elec_dat(i, j, b, :, r)), 'Method', 'fmcd');
                ticker = ticker + 1;
                waitbar(ticker / 192000, h);
            end
        end
    end
end

df = size(dd, 1) - 1;

robust_tstat = mu ./ (sqrt(sig ./ size(dd, 1)));
robust_tstat = mu ./ (sqrt(sig ./ 25));
tdist2T = @(t,v) (1 - betainc(v ./ (v + t .^2), v ./ 2, 0.5));    % 2-tailed t-distribution
for i = 1:6
    for j = 1:5
%        p(:, :, i, j) = tdist2T(robust_tstat(:, :, i, j), 25 - 1);
         p(:, :, i, j) = 1 - tcdf((robust_tstat(:, :, i, j)), 25 - 1); %% change to correct dof
    end
end


%% Perform FDR across time in all electrodes per beta
time_idx = 35:104;
p = p(:, time_idx, :, :);

% correct cluster sizes within electrode
cluster_pval = 0.05;
voxel_pval = 0.001;
clustsizes = [];
for r = 1:6
    for e = 1:5
        clustinfo = bwconncomp(p(:, :, r, e) < voxel_pval);
        clustsizes = [clustsizes; cellfun(@numel, clustinfo.PixelIdxList)'];
        
        clust_threshold(r, e) = prctile(cellfun(@numel, clustinfo.PixelIdxList), 100 - cluster_pval * 100);
    end
end


for r = 1:6
    for e = 1:5
        clustinfo = bwconncomp(p(:, :, r, e) < voxel_pval);
        % identify clusters to remove
        clust_info = cellfun(@numel, clustinfo.PixelIdxList);
        whichclusters2remove = find(clust_info < clust_threshold(r, e));
        
        % remove clusters
        ptemp = p(:, :, r, e) < voxel_pval;
        for i = 1:length(whichclusters2remove)
            ptemp(clustinfo.PixelIdxList{whichclusters2remove(i)}) = 0;
        end
        p_indiv(:, :, r, e) = ptemp;
    end
end
disp('here are some example plots')
figure
contourf(timerstamps(time_idx), freqz, p_indiv(:, :, 1, 3) .* mu(:, time_idx, 1, 3), 'linecolor', 'none');
title('cluster corrected mean b0 at electrode C3');
figure;
contourf(timerstamps(time_idx), freqz, p_indiv(:, :, 1, 2) .* mu(:, time_idx, 1, 2), 'linecolor', 'none');
title('cluster corrected mean b0 at electrode Fcz');
figure;
contourf(timerstamps(time_idx), freqz, p_indiv(:, :, 4, 2) .* mu(:, time_idx, 4, 2), 'linecolor', 'none');
title('cluster corrected mean Error*IL at electrode Fcz');




%% I STOPPEd HERE
%% THIS IS A THRESHOLD OF CLUSTER SIZE USING CLUSTERS FROM ALL ELECTRODES
all_elec_thres = round(prctile(clustsizes, 100 - cluster_pval * 100));

%%%%%% END USING CLUSTER CORRECTION APPROACH %%%%%%


%% You can just correct at particualr frequency bands using FDR
theta_pval = squeeze(mean(p(2:5, :, :, :), 1));

beta_pval = squeeze(mean(p(12:17, :, :, :), 1));

beta_high_pval = squeeze(mean(p(18:28, :, :, :), 1));
error_rate = 0.1;
for r = 1:6
    for e = 1:5
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

for r = 1:6
    for e = 1:5
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


for r = 1:6
    for e = 1:5
        p1 = beta_high_pval(:, r, e);
            p1 = fdr_bh(p1(:)', error_rate);
            beta_high_h1(:, r, e) = p1';
    end
end
figure;
imagesc(timerstamps(time_idx), 1:6, beta_high_h1(:, :, 3)');
ylabel('regressor')
xlabel('time (ms)')
title('betalow(20-30 Hz) activity across regressors in electode C3')


%% WHOLE FDR TxF within electrode
for r = 1:6
    for e = 1:5
        
        p1 = p(:, :, r, e);
        p1 = fdr_bh(p1(:)', 0.15);
        fp_all(:, :, r, e) = reshape(p1, 32, size(p, 2));
        
    end
end

figure;
contourf(timerstamps(time_idx), freqz, fp_all(:, :, 4, 2) .* mu(:, time_idx, 4, 2), 'Linecolor', 'none');
title('example of C3 Error*IL regressor using whole electrode Time x frequency FDR')
%% DIFFERENT APPROACH: set cluster threshold to k=20 and use FDR on all
%%% NOT COMPLETED

%% YOU COULD ALSO FIND THE PEAK WITHIN A TIME AND FREQUENCY WINDOW