function compareICA2(EEG, c1, c2, t1, t2)
% compareICA2({EEG, EEG}, 16, 33, [50, 80], [62, 92])
% c1: IC 1, c2: IC 2
% t1: window for IC 1
% t2: window for IC 2
c1d = 1; c2d = 1;
if c1 < 0
    c1 = abs(c1);
    c1d = -1;
end
if c2 < 0
    c2 = abs(c2);
    c2d = -1;
end

time_win = [-500, 500]; % -500 ms to 500 ms

if numel(EEG) == 1
    figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    %{
    % Does this section remove the EOG channel due to running ICA without removing EOG channel?
    kkkk = size(EEG.data);
    EEG.data = reshape(EEG.data, 64, []);
    EEG.icaact = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :);
    EEG.data = reshape(EEG.data, kkkk);
    EEG.icaact = reshape(EEG.icaact, kkkk - [1, 0, 0]);
    
    % kkkk=size(FF.data);
    % FF.data=reshape(FF.data,64,[]);
    % FF.icaact = (FF.icaweights*FF.icasphere)*FF.data(FF.icachansind,:);
    % FF.data=reshape(FF.data,kkkk);
    % FF.icaact=reshape(FF.icaact,kkkk);
    %}
    
    xlim = dsearchn(EEG.times', time_win');
    subplot(2, 5, 1)
    topoplot(c1d * EEG.icawinv(:, c1), EEG.chanlocs(EEG.icachansind));
    % kkkkk is the labeled class
    [~, kkkkk] = max( EEG.etc.ic_classification.ICLabel.classifications(c1, :) );
    title( {EEG.dipfit.model(c1).areadk, EEG.setname, [num2str(c1), ' ', EEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * EEG.etc.ic_classification.ICLabel.classifications(c1, kkkkk), '%2.1f')]} )
    
    subplot(2, 5, 2)
    topoplot(c1d * EEG.dipfit.model(c1).sourcepot, EEG.chanlocs(EEG.icachansind));
    title('source')
    
    subplot(2, 5, 3)
    topoplot(c1d * EEG.dipfit.model(c1).diffmap, EEG.chanlocs(EEG.icachansind));
    title('error')
    
    subplot(2, 5, 6)
    topoplot(c2d * EEG.icawinv(:, c2), EEG.chanlocs(EEG.icachansind));
    [~, kkkkk] = max( EEG.etc.ic_classification.ICLabel.classifications(c2,:) );
    title( {EEG.dipfit.model(c2).areadk, EEG.setname, [num2str(c2), ' ', EEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * EEG.etc.ic_classification.ICLabel.classifications(c2,kkkkk), '%2.1f')]}, 'Color', 'Red')
    
    
    subplot(2, 5, 7)
    topoplot(c2d * EEG.dipfit.model(c2).sourcepot, EEG.chanlocs(EEG.icachansind));
    title('source', 'Color', 'Red')
    subplot(2, 5, 8)
    topoplot(c2d * EEG.dipfit.model(c2).diffmap, EEG.chanlocs(EEG.icachansind));
    title('error', 'Color', 'Red')
    
    subplot(2, 5, 4)
    plot(EEG.times(xlim(1):xlim(2)), c1d * mean(EEG.icaact(c1, xlim(1):xlim(2), :), 3), 'k')
    hold on
    plot(EEG.times(xlim(1):xlim(2)), c2d * mean(EEG.icaact(c2, xlim(1):xlim(2), :), 3), 'r')
    title('Activation')
    vline([t1(1), t1(2)], {'k--', 'k--'}, {num2str(t1(1)), []})
    vline([t2(1), t2(2)], {'r--', 'r--'}, {num2str(t2(1)), []})
    %alpha(h2,.2)
    set(gca, 'xlim', time_win)
    
    subplot(2, 5, 9)
    plot(EEG.times(xlim(1):xlim(2)), c1d * normalize(mean(EEG.icaact(c1, xlim(1):xlim(2), :), 3)),'k')
    hold on
    plot(EEG.times(xlim(1):xlim(2)), c2d * normalize(mean(EEG.icaact(c2, xlim(1):xlim(2), :), 3)),'r')

    title('NORMALIZED Activation')
    vline([t1(1), t1(2)], {'k--', 'k--'})
    vline([t2(1), t2(2)], {'r--', 'r--'})
    set(gca, 'xlim', time_win)
    
    subplot(2, 5, [5, 10])
    nnn = dsearchn(EEG.times', t1');
    x1 = squeeze(squeeze( mean(EEG.icaact(c1, nnn(1):nnn(2), :), 2) ));
    nnn = dsearchn(EEG.times', t2');
    x2 = squeeze(squeeze( mean(EEG.icaact(c2, nnn(1):nnn(2), :), 2) ));
    plot(x1, x2, '.')
    [rr, pp] = corr(x1, x2, 'Type', 'Spearman');
    title(['r across trials: r=', num2str(rr, '%0.3f'), ' p=', num2str(pp,'%0.3f')])
else
    figure
    EEG1 = EEG{1};
    EEG2 = EEG{2};
    
    %{
    kkkk=size(EEG1.data);
    EEG1.data=reshape(EEG1.data,64,[]);
    EEG1.icaact = (EEG1.icaweights*EEG1.icasphere)*EEG1.data(EEG1.icachansind,:);
    EEG1.data=reshape(EEG1.data,kkkk);
    EEG1.icaact=reshape(EEG1.icaact,kkkk-[1 0 0]);
    
    kkkk=size(EEG2.data);
    EEG2.data=reshape(EEG2.data,64,[]);
    EEG2.icaact = (EEG2.icaweights*EEG2.icasphere)*EEG2.data(EEG2.icachansind,:);
    EEG2.data=reshape(EEG2.data,kkkk);
    EEG2.icaact=reshape(EEG2.icaact,kkkk-[1 0 0]);
    %}
    
    xlim = dsearchn(EEG1.times', time_win');
    subplot(2,7,1)
    topoplot(c1d*EEG1.icawinv(:,c1),EEG1.chanlocs(EEG1.icachansind));
    kkkkk=find(EEG1.etc.ic_classification.ICLabel.classifications(c1,:)==max(EEG1.etc.ic_classification.ICLabel.classifications(c1,:)));
    title({EEG1.dipfit.model(c1).areadk,EEG1.setname,[num2str(c1) ' ' EEG1.etc.ic_classification.ICLabel.classes{kkkkk} ' ' num2str(100*EEG1.etc.ic_classification.ICLabel.classifications(c1,kkkkk),'%2.1f') ]})
    
    subplot(2,7,2)
    topoplot(c1d*EEG1.dipfit.model(c1).sourcepot,EEG1.chanlocs(EEG1.icachansind));
    title('source')
    
    subplot(2,7,3)
    topoplot(c1d*EEG1.dipfit.model(c1).diffmap,EEG1.chanlocs(EEG1.icachansind));
    title('error')
    
    subplot(2,7,8)
    topoplot(c2d*EEG2.icawinv(:,c2),EEG2.chanlocs(EEG2.icachansind));
    kkkkk=find(EEG2.etc.ic_classification.ICLabel.classifications(c2,:)==max(EEG2.etc.ic_classification.ICLabel.classifications(c2,:)));
    title({EEG2.dipfit.model(c2).areadk,EEG2.setname,[ num2str(c2) ' ' EEG2.etc.ic_classification.ICLabel.classes{kkkkk} ' ' num2str(100*EEG2.etc.ic_classification.ICLabel.classifications(c2,kkkkk),'%2.1f') ]},'Color','Red')
    
    
    subplot(2,7,9)
    topoplot(c2d*EEG2.dipfit.model(c2).sourcepot,EEG2.chanlocs(EEG2.icachansind));
    title('source','Color','Red')
    
    subplot(2,7,10)
    topoplot(c2d*EEG2.dipfit.model(c2).diffmap,EEG2.chanlocs(EEG2.icachansind));
    title('error','Color','Red')
    
    
    subplot(2,7,4)
    gg = EEG1.dipfit.model(c1).datapot .* EEG1.data(EEG1.icachansind, :, :);
    plot(EEG1.times(xlim(1):xlim(2)), mean(c1d * squeeze(gg(c1, xlim(1):xlim(2), :)), 2), 'k-', 'LineWidth', 0.5)
    hold on
    gg = EEG2.dipfit.model(c2).datapot.*EEG2.data(EEG2.icachansind,:,:);
%    plot(EEG2.times(xlim(1):xlim(2)),mean(c2d*squeeze(gg(c2,xlim(1):xlim(2),:)),2),'r-','LineWidth',0.5)

    plot(EEG1.times(xlim(1):xlim(2)), c1d * mean(EEG1.icaact(c1, xlim(1):xlim(2), :), 3), 'k-', 'LineWidth', 2)
    %gg=EEG1.dipfit.model(c1).diffmap.*EEG1.data;
    %plot(EEG1.times(xlim(1):xlim(2)),mean(c1d*squeeze(gg(c1,xlim(1):xlim(2),:)),2),'k:','LineWidth',1)
    
    plot(EEG2.times(xlim(1):xlim(2)), c2d * mean(EEG2.icaact(c2, xlim(1):xlim(2), :), 3), 'r-', 'LineWidth', 2)
    %gg=EEG2.dipfit.model(c2).diffmap.*EEG2.data;
    %plot(EEG2.times(xlim(1):xlim(2)),mean(c2d*squeeze(gg(c2,xlim(1):xlim(2),:)),2),'r:','LineWidth',1)
    title('Activation')
%     vline([t1(1), t1(2)], {'k--', 'k--'})
%     vline([t2(1), t2(2)], {'r--', 'r--'})
    %alpha(h2,.2)
    set(gca,'xlim',time_win)
    
    subplot(2,7,11)
    plot(EEG1.times(xlim(1):xlim(2)),c1d*normalize(mean(EEG1.icaact(c1,xlim(1):xlim(2),:),3)),'k')
    hold on
    plot(EEG2.times(xlim(1):xlim(2)),c2d*normalize(mean(EEG2.icaact(c2,xlim(1):xlim(2),:),3)),'r')
    title('NORMALIZED Activation')
    vline([t1(1), t1(2)], {'k--', 'k--'})
    vline([t2(1), t2(2)], {'r--', 'r--'})
    set(gca,'xlim',time_win)
    
    nnn=dsearchn(EEG1.times',t1');
    x1=squeeze(squeeze(mean(EEG1.icaact(c1,nnn(1):nnn(2),:),2)));
    nnn=dsearchn(EEG2.times',t2');
    x2=squeeze(squeeze(mean(EEG2.icaact(c2,nnn(1):nnn(2),:),2)));
    
    hh=subplot(2,7,5);
    imagesc(EEG1.times(xlim(1):xlim(2)),1:size(EEG1.data,3),c1d*squeeze(EEG1.icaact(c1,xlim(1):xlim(2),:))')
    dRange = 3*median(median(abs(squeeze(EEG1.icaact(c1,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])

    hh=subplot(2,7,6);
    gg=EEG1.dipfit.model(c1).sourcepot.*EEG1.data(EEG1.icachansind,:,:);
    imagesc(EEG1.times(xlim(1):xlim(2)),1:size(EEG1.data,3),c1d*squeeze(gg(c1,xlim(1):xlim(2),:))')
    dRange = 3*median(median(abs(squeeze(gg(c1,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])

    hh=subplot(2,7,7);
    
    TMP = EEG1;
    %{
    kkkk=size(TMP.data);
    TMP.data=reshape(TMP.data,64,[]);
    TMP.icaweights(:,c1) = TMP.dipfit.model(c1).sourcepot;
    TMP.icaact = (TMP.icaweights*TMP.icasphere)*TMP.data(TMP.icachansind,:);
    TMP.data=reshape(TMP.data,kkkk);
    TMP.icaact=reshape(TMP.icaact,kkkk-[1 0 0]);
    %}
    
%    imagesc(EEG1.times(xlim(1):xlim(2)),1:size(EEG1.data,3),c1d*squeeze(gg(c1,xlim(1):xlim(2),:))')
imagesc(EEG1.times(xlim(1):xlim(2)),1:size(EEG1.data,3),c1d*squeeze(EEG1.icaact(c1,xlim(1):xlim(2),:))'-c1d*squeeze(TMP.icaact(c1,xlim(1):xlim(2),:))')
%    dRange = 3*median(median(abs(squeeze(gg(c1,xlim(1):xlim(2),:))')));
    dRange = 3*median(median(abs(c1d*squeeze(EEG1.icaact(c1,xlim(1):xlim(2),:))'-c1d*squeeze(gg(c1,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])
    
    hh=subplot(2,7,12);
    imagesc(EEG2.times(xlim(1):xlim(2)),1:size(EEG2.data,3),c2d*squeeze(EEG2.icaact(c2,xlim(1):xlim(2),:))')
    dRange = 3*median(median(abs(squeeze(EEG2.icaact(c2,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])
    
    hh=subplot(2,7,13);
    gg=EEG2.dipfit.model(c2).sourcepot .* EEG2.data(EEG2.icachansind,:,:);
    imagesc(EEG2.times(xlim(1):xlim(2)),1:size(EEG2.data,3),c2d*squeeze(gg(c2,xlim(1):xlim(2),:))')
    dRange = 3*median(median(abs(squeeze(gg(c2,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])

    hh=subplot(2,7,14);
    gg=EEG2.dipfit.model(c2).diffmap.*EEG2.data(EEG2.icachansind,:,:);
    imagesc(EEG2.times(xlim(1):xlim(2)),1:size(EEG2.data,3),c2d*squeeze(gg(c2,xlim(1):xlim(2),:))')
    dRange = 3*median(median(abs(squeeze(gg(c2,xlim(1):xlim(2),:))')));
    set(hh,'clim',[-dRange dRange])
end
