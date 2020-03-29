close all; clearvars; clc

%% load aligned data
summarypath = uigetdir;
summaryfile = fullfile(summarypath, 'behavior_resultSummary.mat');
summary = load(summaryfile);

%% alinged with lift onset
All_info_trial = summary.All_info_trial;
nsub = height(summary.All_info_trial);
for sub_i = 1:nsub
    sub_id = summary.All_info_trial.subID{sub_i, 1};
    data = summary.All_info_trial.data{sub_i, 1};
    dt = summary.All_info_trial.dt(sub_i, 1); % in ms
    data_lft_onset_aligned = summary.All_info_trial.data{sub_i, 1}.onset_aligned;
    [~, lft_onset] = min(abs(data_lft_onset_aligned{1, 1}.time)); % lift onset is at 0 s
    nep = height(data);
    
    %% find the distance between the dCOPy, dFy, Fgrip to the ideal surface
    dCOPz = -29.6 * 2; % from VF to TH
    dist = nan(nep, 1);
    for ep = 1:nep
        if strcmpi(data.side(ep, 1), 'R')
            Mcom = 395;
        else
            Mcom = -395;
        end
        
        dCOPy = data_lft_onset_aligned{ep, 1}.dCOPy(lft_onset);
        dFy = data_lft_onset_aligned{ep, 1}.dFy(lft_onset);
        Fgrip = data_lft_onset_aligned{ep, 1}.GripF(lft_onset);
        
        dist(ep, 1) = dist2surface(Mcom, dCOPz, dCOPy, dFy, Fgrip);
    end
    
    dist = array2table(dist, 'VariableNames', {'dist2ideal'});
    All_info_trial.data{sub_i, 1} = horzcat(All_info_trial.data{sub_i, 1}, dist);
end
save(fullfile(summarypath, 'behavior_resultSummary_w_dist2ideal.mat'), 'All_info_trial');
