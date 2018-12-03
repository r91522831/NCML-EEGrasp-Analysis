function [eegout] = checkEvent(eegin, nbevent)
%checkEvent Summary of this function goes here
%   Detailed explanation goes here
eegout = eegin;

start_count = false;
% s9 is the trigger for trial beginning
% s129 is the trigger for trial end
event_check = strcmp({eegin.event.type}, 's9') + 2 * strcmp({eegin.event.type}, 's129');
keep_event = true(size(event_check));
for i = 1:length(event_check)
    if event_check(i) == 1
        start_count = true;
        trial_on = i;
        event_count = 1;
        continue;
    end
    
    if event_check(i) == 2
        if i - trial_on < nbevent - 1
            for k = trial_on:i
                keep_event(k) = false;
            end
        end
        
        event_count = 0;
        start_count = false;
        continue;
    end
    
    if start_count
        event_count = event_count + 1;
    end
end

eegout.event = eegout.event(keep_event);

end

