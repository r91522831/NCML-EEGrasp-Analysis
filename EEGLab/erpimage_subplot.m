function erpimage_subplot(data, freqName, chanName, col)
% erpimage_subplot Summary of this function goes here
%   subplot for erp images based on frequency and channel
%   Detailed explanation goes here

chanID = find(strcmpi({data.chanlocs.labels}, chanName));

subplot(3, 3, 6 + col);
condType = 'IL';
tmp_trial_id = strcmp({data.epoch.condType}, condType);
n_trial = length(find(tmp_trial_id));
[~, ~, ~, ~, axh] = erpimage(squeeze(data.data(chanID, :, tmp_trial_id)), [], data.times, [], floor(n_trial / 3), [], 'erp', 'on');
axh{1}.YLabel.String = condType;

subplot(3, 3, 3 + col);
condType = 'TR';
tmp_trial_id = strcmp({data.epoch.condType}, condType);
n_trial = length(find(tmp_trial_id));
[~, ~, ~, ~, axh] = erpimage(squeeze(data.data(chanID, :, tmp_trial_id)), [], data.times,[], floor(n_trial / 3), [], 'erp', 'on');
axh{1}.XTickLabels = [];
axh{1}.XLabel = [];
axh{1}.YLabel.String = condType;
axh{3}.XTickLabels = [];
axh{3}.XLabel = [];

subplot(3, 3, col);
condType = 'PT';
tmp_trial_id = strcmp({data.epoch.condType}, condType);
n_trial = length(find(tmp_trial_id));
[~, ~, ~, ~, axh] = erpimage(squeeze(data.data(chanID, :, tmp_trial_id)), [], data.times, [], floor(n_trial / 3), [], 'erp', 'on');
axh{1}.XTickLabels = [];
axh{1}.XLabel = [];
axh{1}.XAxisLocation = 'top';
axh{1}.XLabel.String = [freqName, ' ', data.chanlocs(chanID).labels];
axh{1}.YLabel.String = condType;
axh{3}.XTickLabels = [];
axh{3}.XLabel = [];
end

