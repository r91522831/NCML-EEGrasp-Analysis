figure
sub = {'S002 Aug 2018', 'S004 Aug 2018', 'S005 Aug 2018', 'S006 Aug 2018', ...
       'S007 Feb 2019', ...
       'S009 Aug 2018', ...
       'S010 Feb 2019', 'S011 Feb 2019', 'S012 Feb 2019', 'S013 Feb 2019', ...
       'S014 Sep 2019', 'S015 Sep 2019', ...
       'S016 Feb 2019', 'S017 Feb 2019', ...
       'S018 Sep 2019', 'S019 Sep 2019', 'S021 Sep 2019', 'S022 Sep 2019', 'S023 Sep 2019', 'S024 Sep 2019', 'S025 Sep 2019', 'S026 Sep 2019', 'S027 Sep 2019'};
for i = 1:23
    subplot(6, 4, i)
    plot(all_sub_mx_pRoll(:, 1, i), all_sub_mx_pRoll(:, 2, i), '.')
    xlabel(sub{i})
    axis([-600, 600, -25, 25])
end


%%
non_zero = [true(1); info_time_trigger{i, 1}(2:end) ~= 0];
figure
subplot(3, 1, 1)
plot(info_time_trigger{i, 1}(non_zero)./1000, resultantF{i, 1}{non_zero, 'fy'});
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 6})./1000, ':r');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000, ':b');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 2})./1000, '--k');
subplot(3, 1, 2)
plot(info_time_trigger{i, 1}(non_zero)./1000, obj_height{i, 1}(non_zero));
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 6})./1000, ':r');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000, ':b');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 2})./1000, '--k');
subplot(3, 1, 3)
for j = 1:length(angTilt2R{i, 1})
    if angTilt2R{i, 1}(j) > 90
        angTilt2R{i, 1}(j) = angTilt2R{i, 1}(j) - 180;
    elseif angTilt2R{i, 1}(j) < -90
        angTilt2R{i, 1}(j) = angTilt2R{i, 1}(j) + 180;
    end
end
plot(info_time_trigger{i, 1}(non_zero)./1000, angTilt2R{i, 1}(non_zero));
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 6})./1000, ':r');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000, ':b');
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 2})./1000, '--k');