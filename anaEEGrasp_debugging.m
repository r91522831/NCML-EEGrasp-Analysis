figure
subplot(3, 1, 1)
plot(info_time_trigger{i, 1}./1000, resultantF{i, 1}.mx)
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000)
subplot(3, 1, 2)
plot(info_time_trigger{i, 1}./1000, obj_height{i, 1})
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000)
subplot(3, 1, 3)
plot(info_time_trigger{i, 1}./1000, angTilt2R{i, 1})
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000)




figure
subplot(3, 1, 1)
plot(info_time_trigger{i, 1}(1:end-1)./1000, resultantF{i, 1}{1:end-1, 'mx'});
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000);
subplot(3, 1, 2)
plot(info_time_trigger{i, 1}(1:end-1)./1000, obj_height{i, 1}(1:end-1));
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000);
subplot(3, 1, 3)
for j = 1:length(angTilt2R{i, 1})
    if angTilt2R{i, 1}(j) > 90
        angTilt2R{i, 1}(j) = angTilt2R{i, 1}(j) - 180;
    elseif angTilt2R{i, 1}(j) < 90
        angTilt2R{i, 1}(j) = angTilt2R{i, 1}(j) + 180;
    end
end
plot(info_time_trigger{i, 1}(1:end-1)./1000, angTilt2R{i, 1}(1:end-1));
vline(info_time_trigger{i, 1}(ind_lft_onset{i, 3})./1000);