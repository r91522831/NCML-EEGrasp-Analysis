clear;
%% Insert behavior onset into EEG evenets for SPM
[filename, pathname, ~] = uigetfile;
% load('/Users/yenhsunw/Dropbox (ASU)/NCML-EEGrasp/EEG/SPM/S009_spm_raw.mat')
load(fullfile(pathname, filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get experiment conditions
sub_dir = uigetdir('', 'select subject folder for behvior matlab raw data.');
file_list = dir(fullfile(sub_dir, '*.csv'));
tmp = char({file_list.name});
cond = tmp(:, [11:12, 16:18]);

for i = 1:length(D.trials)
    D.trials(i).label = cond(i, :);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(pathname, filename), 'D');
