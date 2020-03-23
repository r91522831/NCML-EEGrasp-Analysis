close all; clear; clc

res = load('/Users/yenhsunw/Dropbox (ASU)/BIDS_format/NCML-EEGraspU_pilot/behavior_resultSummary.mat');

%%
sub_i = 3;
tbl = res.All_info_trial.data{sub_i, 1}(:, {'dCOPy', 'dFy', 'gripF', 'pRoll'});
mdl = fitlm(tbl);
disp(mdl)
