close all; clear; clc;

% LIMO level 1 analysis
Y = EEG.data;
[nb_ch, nb_time, nb_epoc] = size(Y);

X = double([ strcmp([EEG.epoch.cond]', {'IL'}), ...
             strcmp([EEG.epoch.cond]', {'TR'}), ...
             strcmp([EEG.epoch.cond]', {'PT1'}),...
             strcmp([EEG.epoch.cond]', {'PT2'}),...
             strcmp([EEG.epoch.cond]', {'PT3'}) ]);
nb_cond = size(X, 2);
         
B = nan(nb_ch, nb_cond);
for i = 1:nb_ch
    B(i, :) = diag(pinv(X' * X) * X' * squeeze(Y(i, :, :))');
end
