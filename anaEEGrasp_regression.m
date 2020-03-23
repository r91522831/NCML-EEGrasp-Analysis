close all; clear; clc

res = load('/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGraspU_pilot/behavior_resultSummary.mat');

%%
for sub_i = 1:3
    tbl = res.All_info_trial.data{sub_i, 1}(:, {'dCOPy', 'dFy', 'gripF', 'pRoll'});
    mdl = fitlm(tbl);
    disp(res.All_info_trial.subID{sub_i, 1})
    disp(mdl)
    disp('===============================================================')
end


