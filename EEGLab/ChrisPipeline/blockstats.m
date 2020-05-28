function [avgBlock, seBlock] = blockstats(data, epblock)
%blockstats Summary of this function goes here
%   Detailed explanation goes here
nepb = length(epblock);
[nIC, ntwin_eeg, nep] = size(data);
% nIC x ntwin_eeg x nepb
avgBlock = nan(nIC, ntwin_eeg, nepb);
seBlock = avgBlock;
for i_epb = 1:nepb
    avgBlock(:, :, i_epb) = mean(data(:, :, epblock{1, i_epb}), 3);
    seBlock(:, :, i_epb) = std(data(:, :, epblock{1, i_epb}), 0, 3) ./ sqrt(length(epblock{1, i_epb}));
end

% for all epoch
avgBlock(:, :, nepb + 1) = mean(data(:, :, :), 3);
seBlock(:, :, nepb + 1) = std(data(:, :, :), 0, 3) ./ sqrt(nep);
end

