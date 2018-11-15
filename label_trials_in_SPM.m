clear;
%% Label trial conditions for SPM
[filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009_spm_raw.mat')
load(fullfile(pathname, filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get experiment conditions
sub_dir = uigetdir('', 'select subject folder for behvior matlab raw data.');
file_list = dir(fullfile(sub_dir, '*.csv'));
tmp = char({file_list.name});

cond = tmp(:, [11:12, 16:18]);
session = [tmp(:, 11:12), num2str(zeros(length(tmp), 1)), tmp(:, 13:14)];

for i = 1:length(D.trials)
    if strcmp(cond(i, 1), 'I')
        D.trials(i).label = session(i, :);
        D.trials(i).tag = cond(i, :);
    else
        D.trials(i).label = cond(i, :);
        D.trials(i).tag = session(i, :);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(pathname, filename), 'D');
