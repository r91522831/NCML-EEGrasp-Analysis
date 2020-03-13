
figure
% % % info_time_trigger{ep, 1}();
hold on
for ep = 1:length(obj_height)
    
    [~, lft_ind] = min(abs(info_onset_time(ep, 1) * 1000 - info_time_trigger{ep, 1}));
    tmp_ind = (lft_ind + (-500:500))';
    plot(0.001 * (-500:500), obj_height{ep, 1}(tmp_ind, 1))
end

vline(0, '--r', 'onset')
ylabel('obj height (mm)')
xlabel('time (s)')
set(gca, 'FontSize', 18)