function plot_preliminary(band, e1, e2, f1, f2, t1, t2, electrodes_name, freqz, timerstamps, mu, err_range)
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
err = [max(err_range(1, selected_t)); min(err_range(2, selected_t))];
ersp_est_IL = mean(selected_beta0(:)) + mean(selected_beta3(:)) * err; % in dB
ersp_est_TR = mean(selected_beta1(:)) + mean(selected_beta4(:)) * err; % in dB
ersp_est_PT = mean(selected_beta2(:)) + mean(selected_beta5(:)) * err; % in dB

% dB = 10 log(1 + X); X = 0.01 -> 1%
% db2pow: y dB = 10 log (y)
plot(err, 100 * (db2pow(ersp_est_IL) - 1), '-ob', err, 100 * (db2pow(ersp_est_TR) - 1), '-sr', err, 100 * (db2pow(ersp_est_PT) - 1), '--dr');
ylabel('estERSP (%)')
xlabel({'err ({\circ})', [band, ': ', e1, ' and ', e2, ', freq: ', num2str(f1), ' to ', num2str(f2), ' Hz, time: ', num2str(t1), ' to ', num2str(t2), ' ms']})
end

