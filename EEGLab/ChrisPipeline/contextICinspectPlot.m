function fig = contextICinspectPlot(EEG, data, context, time, blockavg, blockse, comp, plotvar, powfft, freq, powblockavg, powblockse)
% context = {'IL_1', 'IL_{2-19}', 'TR_1', 'TR_{2-19}', 'PT1_1', 'PT1_{2-19}', 'PT2_1', 'PT2_{2-19}', 'PT3_1', 'PT3_{2-19}', 'all'};
% % % 'LineStyle', 'none', 'Color'
linspec = {'-', 'none', '-', 'none', '-', 'none', [], 'none', [], 'none', 'none'};
lincolor = {'r', 'r', 'b', 'b', 'k', 'k', [], 'r', [], 'b', 'k'};
nepb = length(context);
dt = diff(0.001 * time);
cutoff = 15; % Hz

fig = figure('DefaultAxesFontSize', 18, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
subplot(3, 3, 1)
tmp = [ data(:, :, context{1, 1}), data(:, :, context{1, 2}(1, 1)), ...
        data(:, :, context{1, 3}), data(:, :, context{1, 5}) ];
dRange = nanmax(abs(tmp), [], 'all');
hold on
plot(0.001 * time, filtmat_class( dt(1), cutoff, data(:, :, context{1, 1})' )); % IL1
plot(0.001 * time, filtmat_class( dt(1), cutoff, data(:, :, context{1, 2}(1, 1))' )); % IL2
plot(0.001 * time, filtmat_class( dt(1), cutoff, data(:, :, context{1, 3})' )); % TR1
plot(0.001 * time, filtmat_class( dt(1), cutoff, data(:, :, context{1, 5})' )); % PT1 1
hold off
ylim([-dRange, dRange]);
xlim(0.001 * [time(1), time(end)]);
vline(0, '--k');
if strcmpi(plotvar, 'icaact'), ylabel('\muV'); else, ylabel('10log10(\muV^{2})'); end

legend({'IL_1', 'IL_2', 'TR_1', 'PT1_1'}, 'Location', 'best')
title(plotvar)

subplot(3, 3, 4)
hold on
for i_epb = 2:2:6
    h(i_epb) = shadedErrorBar(0.001 * time, filtmat_class( dt(1), cutoff, blockavg(:, :, i_epb)' ), blockse(:, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
end
hold off
% % %     ylim([-dRange, dRange]);
xlim(0.001 * [time(1), time(end)]);
vline(0, '--k');
if strcmpi(plotvar, 'icaact'), ylabel('\muV'); else, ylabel('10log10(\muV^{2})'); end
tmp_shade = [h(2:2:6).patch];
legend(tmp_shade, {'IL_{2~19}', 'TR_{2~19}', 'PT1_{2~19}'}, 'Location', 'best')

subplot(3, 3, 7)
ep_all = nepb + 1;
h(ep_all) = shadedErrorBar(0.001 * time, filtmat_class( dt(1), cutoff, blockavg(:, :, ep_all)' ), blockse(:, :, ep_all), {'LineStyle', linspec{ep_all}, 'Color', lincolor{ep_all}}, 1);
hold on
for i_epb = 8:2:10
    h(i_epb) = shadedErrorBar(0.001 * time, filtmat_class( dt(1), cutoff, blockavg(:, :, i_epb)' ), blockse(:, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
end
hold off
% % %     ylim([-dRange, dRange]);
xlim(0.001 * [time(1), time(end)]);
vline(0, '--k');
if strcmpi(plotvar, 'icaact'), ylabel('\muV'); else, ylabel('10log10(\muV^{2})'); end
xlabel('time (s)')
tmp_shade = [h(ep_all).patch, h(8).patch, h(10).patch];
legend(tmp_shade, {'all', 'PT2_{2~19}', 'PT3_{2~19}'}, 'Location', 'best')


subplot(3, 3, 2)
tmp = [ powfft(:, :, context{1, 1}), powfft(:, :, context{1, 2}(1, 1)), ...
        powfft(:, :, context{1, 3}), powfft(:, :, context{1, 5}) ];
dRange = nanmax(abs(tmp), [], 'all');
hold on
plot(freq, powfft(:, :, context{1, 1})); % IL1
plot(freq, powfft(:, :, context{1, 2}(1, 1))); % IL2
plot(freq, powfft(:, :, context{1, 3})); % TR1
plot(freq, powfft(:, :, context{1, 5})); % PT1 1
hold off
ylim([0, dRange]);
xlim([freq(1), freq(end)]);
vline(0, '--k');
ylabel('\muV^2');

legend({'IL_1', 'IL_2', 'TR_1', 'PT1_1'}, 'Location', 'best')
title('power spectrum')

subplot(3, 3, 5)
hold on
for i_epb = 2:2:6
    h(i_epb) = shadedErrorBar(freq, powblockavg(:, :, i_epb), powblockse(:, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
end
hold off
% % %     ylim([-dRange, dRange]);
xlim([freq(1), freq(end)]);
vline(0, '--k');
ylabel('\muV^2');
tmp_shade = [h(2:2:6).patch];
legend(tmp_shade, {'IL_{2~19}', 'TR_{2~19}', 'PT1_{2~19}'}, 'Location', 'best')

subplot(3, 3, 8)
ep_all = nepb + 1;
h(ep_all) = shadedErrorBar(freq, powblockavg(:, :, ep_all), powblockse(:, :, ep_all), {'LineStyle', linspec{ep_all}, 'Color', lincolor{ep_all}}, 1);
hold on
for i_epb = 8:2:10
    h(i_epb) = shadedErrorBar(freq, powblockavg(:, :, i_epb), powblockse(:, :, i_epb), {'LineStyle', linspec{i_epb}, 'Color', lincolor{i_epb}}, 1);
end
hold off
% % %     ylim([-dRange, dRange]);
xlim([freq(1), freq(end)]);
vline(0, '--k');
ylabel('\muV^2');
xlabel('freq (Hz)')
tmp_shade = [h(ep_all).patch, h(8).patch, h(10).patch];
legend(tmp_shade, {'all', 'PT2_{2~19}', 'PT3_{2~19}'}, 'Location', 'best')












subplot(3, 3, 3:3:9)
topoplot(EEG.icawinv(:, comp), EEG.chanlocs(EEG.icachansind));
[~, kkkkk] = max( EEG.etc.ic_classification.ICLabel.classifications(comp, :) );
tt = { [ num2str(comp), ' ', EEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * EEG.etc.ic_classification.ICLabel.classifications(comp, kkkkk), '%2.1f') ], ...
    EEG.dipfit.model(comp).areadk, ...
    EEG.setname(1:6) };
title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0])


end