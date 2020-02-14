function EEG = putback_nonEEG_v1(newEEG, originEEG, data_mask)
%putback_nonEEG Summary of this function goes here
%   Detailed explanation goes here
% put non-EEG channgel back to the EEG structure
EEG = newEEG;

if nargin < 3
    data_mask = newEEG.etc.clean_sample_mask;
end

switch length(size(originEEG.data))
    case 2  % raw data length
        EEG.data = [newEEG.data; originEEG.data(~strcmp({originEEG.chanlocs.type}, 'EEG'), data_mask)];
    case 3  % data after epoch
        data2d = reshape(originEEG.data, size(originEEG.data, 1), size(originEEG.data, 2) * size(originEEG.data, 3));
        new_data2d = reshape(newEEG.data, size(newEEG.data, 1), size(newEEG.data, 2) * size(newEEG.data, 3));
        data_attached = [new_data2d; data2d(~strcmp({originEEG.chanlocs.type}, 'EEG'), :)];
        EEG.data = reshape(data_attached, size(data_attached, 1), size(newEEG.data, 2), size(newEEG.data, 3));
    otherwise
end

EEG.nbchan = size(EEG.data, 1);
if ~isempty(EEG.chanlocs)
	EEG.chanlocs = [newEEG.chanlocs, originEEG.chanlocs(~strcmp({originEEG.chanlocs.type}, 'EEG'))];
end
EEG = eeg_checkset( EEG );
end

