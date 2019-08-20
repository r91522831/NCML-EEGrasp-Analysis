function plot_preliminary(band, e1, e2, f1, f2, t1, t2, electrodes_name, freqz, timerstamps, mu, sig, err_range)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
selected_e = strcmp(e1, electrodes_name) | strcmp(e2, electrodes_name);
selected_f = freqz >= f1 & freqz <= f2;
selected_t = timerstamps >= t1 & timerstamps <= t2;
selected_beta0 = mu(selected_f, selected_t, 1, selected_e);
selected_beta3 = mu(selected_f, selected_t, 4, selected_e);
selected_beta1 = mu(selected_f, selected_t, 2, selected_e);
selected_beta4 = mu(selected_f, selected_t, 5, selected_e);
selected_beta2 = mu(selected_f, selected_t, 3, selected_e);
selected_beta5 = mu(selected_f, selected_t, 6, selected_e);
selected_beta0_sig = sig(selected_f, selected_t, 1, selected_e);
selected_beta3_sig = sig(selected_f, selected_t, 4, selected_e);
selected_beta1_sig = sig(selected_f, selected_t, 2, selected_e);
selected_beta4_sig = sig(selected_f, selected_t, 5, selected_e);
selected_beta2_sig = sig(selected_f, selected_t, 3, selected_e);
selected_beta5_sig = sig(selected_f, selected_t, 6, selected_e);

err = [max(err_range(1, selected_t)); min(err_range(2, selected_t))];

% Ref:
% Mean and Variance of Linear Combinations (https://newonlinecourses.science.psu.edu/stat414/node/166/)
% Mean and Variance of Sample Mean (https://newonlinecourses.science.psu.edu/stat414/node/167/)
ersp_est_IL = mean(selected_beta0(:)) + mean(selected_beta3(:)) * err; % in dB
ersp_est_TR = mean(selected_beta1(:)) + mean(selected_beta4(:)) * err; % in dB
ersp_est_PT = mean(selected_beta2(:)) + mean(selected_beta5(:)) * err; % in dB

ersp_est_IL_std = sqrt( mean( selected_beta0_sig(:) ) + mean( selected_beta3_sig(:) ) * (err.^2) );  % in dB
ersp_est_TR_std = sqrt( mean( selected_beta1_sig(:) ) + mean( selected_beta4_sig(:) ) * (err.^2) );  % in dB
ersp_est_PT_std = sqrt( mean( selected_beta2_sig(:) ) + mean( selected_beta5_sig(:) ) * (err.^2) );  % in dB

% dB = 10 log(1 + X); X = 0.01 -> 1%
% db2pow: y dB = 10 log (y)
hold on
errorbar(rad2deg(err) - 0.2, 100 * (db2pow(ersp_est_IL) - 1), 100 * (db2pow(ersp_est_IL_std) - 1), '-ob');
errorbar(rad2deg(err)      , 100 * (db2pow(ersp_est_TR) - 1), 100 * (db2pow(ersp_est_TR_std) - 1), '-sr');
errorbar(rad2deg(err) + 0.2, 100 * (db2pow(ersp_est_PT) - 1), 100 * (db2pow(ersp_est_PT_std) - 1), '--dr');
hold off

ylabel('estERSP (%)')
xlabel({'err ({\circ})', [band, ': ', e1, ' and ', e2, ', freq: ', num2str(f1), ' to ', num2str(f2), ' Hz, time: ', num2str(t1), ' to ', num2str(t2), ' ms']})
end

