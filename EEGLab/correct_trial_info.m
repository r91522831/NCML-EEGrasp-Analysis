%%
orginALLEEG = ALLEEG;

for i = 1:length(ALLEEG)
    for j = 1:length(ALLEEG(i).event)
        ALLEEG(i).event(j).trialID = orginALLEEG(i).event(j).condID;
        switch ALLEEG(i).event(j).cond
            case 'IL'
                ALLEEG(i).event(j).condID = '01';
            otherwise
                ALLEEG(i).event(j).condID = num2str(floor(str2double(orginALLEEG(i).event(j).trialID) / 2), '%.2d');
        end
    end
end