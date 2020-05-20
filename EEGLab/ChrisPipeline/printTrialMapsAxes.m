function handlesICA = printTrialMapsAxes(EEG, myxlim, dataType, freq, myLayout, comp, basewin)
% Example usage
% prettyERPs(ERP, [1 2], 3, {[0 200 0], [0 0 255]})
% Use the data in ERP to plot bins 1 and 2
% The standard error of the difference is 3.
% Plot Bin 1 in green   {[0 200 0],
% and  Bin 2 in blue     [0 0 255]}
%
% Note: This works best if
%       (1) each subject has a 2 minus 1 difference wave and
%       (2) ALLERP(4) is the average of many subjects and
%       (3) You've set the grand average function to measure standard error
%
%c = findobj('tag','prettyFigure');
%if ~isempty(c), close(c), end
%c = findobj('tag','prettyFigureTopo');
%if ~isempty(c), close(c), end
h = figure;% clf
set(h, 'Position', [0, 100, 1964, 1555] * .8, 'tag', 'prettyFigure')

if ~isempty(freq)
    xlim = dsearchn(EEG.icatf.tf_times', myxlim');
    X = EEG.icatf.tf_times(xlim(1):xlim(2));
else
    xlim = dsearchn(EEG.times', myxlim');
    X = EEG.times(xlim(1):xlim(2));
end

%if strcmpi(dataType,'ICA')
%else
data = EEG.data(:, xlim(1):xlim(2), :);
DATA = EEG.data(:, xlim(1):xlim(2), :);
%end

%X = EEG.times;
%absMax = max(max(max(data(1:EEG.nbchan, :, :)))) + 1;
%absMin = min(min(min(data(1:EEG.nbchan, :, :)))) - 1;
if numel(EEG.chanlocs) < 40
    figureSize = 3 / numel(EEG.chanlocs);
    positionScale = .8;
    EEG.chanlocs = EEG.chanlocs(1:30);
else
    figureSize = 7 / numel(EEG.chanlocs);
    positionScale = .65;
end
%figureSize = 4.2 / numel(EEG.chanlocs);
figureSize = figureSize * .9;
positionScaleX = .70;
positionScaleY = .80;
imscProp = .83;
padding = .04;
pad2 = .1;
%
nComps = size(data, 1);
if strcmpi(dataType, 'ICA')
% % %     nComps = size(EEG.icawinv, 2);
    nComps = length(comp);
end
if nComps > 80, nComps=80; end
ncols = floor(nComps / 5);
if mod(ncols,2) == 1, ncols = ncols + 1; end
nrows = 2 * (floor(2 * nComps/ncols) + 1);

%ncols=ncols*2;
nrow = myLayout(1);
ncol = myLayout(2);
if strcmpi(dataType,'ICA')
    plotIdx = -1; pFix=0;
    %for idx=1:nComps
% % %     idx = comp-1;
    fidx = 0;
    for irow = 1:nrow % most likely 5 rows per page
        for icol = 1:ncol
% % %             if idx < nComps
% % %                 idx=idx+1;
            if fidx < nComps
                fidx = fidx + 1;
                idx = comp(fidx);
                
                % ICA activity stacked trial by trial from top to bottom (trial 1 to trial end)
                % [left bottom width height]
                handlesICA(fidx, 1) = axes('position', ...
                                          [ padding / 2 + (icol - 1) / (ncol + pad2), ...
                                            padding / 2 + (nrow - irow) / (nrow + pad2), ...
                                            .95 / (ncol * 2 + pad2) * (1 - padding), ...
                                            .95 / (nrow + pad2) * imscProp * (1 - padding) ]);
                if ~isempty(freq)
                    %{
                    [wt, f] = cwt(reshape(data(idx, :, :), 1, []), EEG.srate);
                    freqidx = dsearchn(f, freq');
                    data = reshape(mean(wt(freqidx(2):freqidx(1), :)), EEG.srate * 4 - 1, []);
                    data = abs(data);
                    baseidx = dsearchn(EEG.times',[-400 0]');
                    %            blc = mean(data(blcidx(1):blcidx(2),:));
                    newData=10*log10( bsxfun(@rdivide, data(:,:), mean(data(baseidx(1):baseidx(2),:),1)));
                    %}
                    freqidx = dsearchn(EEG.icatf.tf_freqs', freq');
                    data = EEG.icatf.tf_ersp{idx}(freqidx(1):freqidx(2), :, :);
                    
                    baseidx = dsearchn(EEG.icatf.tf_times', basewin');
                    basedata = nanmean(data(:, baseidx(1):baseidx(2), :), 2); % baseline for each frequency
                    
                    newData = squeeze(mean(data - basedata, 1)); % average freq band after baseline division for each frequency in '10 * log10( 1^{2}/Hz )'
                    erpY = nanmax(abs(nanmean(newData, 2)));
                    
                    dRange = nanmean(prctile(abs(newData),95)); % dRange=10;
                    imagesc(X, 1:EEG.trials, newData')
                    set(handlesICA(fidx,1),'clim',[-dRange dRange])
                    %title(num2str(dRange))
                else
                    data = EEG.icaact(:, xlim(1):xlim(2), :);
                    erpY = nanmax(abs(nanmean(squeeze(data(idx, :, :)), 2)));
                
                    dRange = nanmean(prctile(abs(data(idx, :, :)), 95));
                    imagesc(X, 1:EEG.trials, squeeze(data(idx, :, :))')
                    set(handlesICA(fidx,1),'clim',[-dRange dRange])
                    if erpY > 1
                        colormap summer
                    else
                        colormap default
                    end
                    %title(num2str(dRange))
                end
                %        set(gca,'XTickLabel',X)
                axis off
                handlesICA(fidx,2)=axes('position',...
                    [padding/2+(icol-1)/(ncol+pad2)
                    padding/2+(nrow-irow)/(nrow+pad2)+.95/(nrow+pad2) * imscProp * (1-padding)
                    .95/(ncol*2+pad2) * (1-padding)
                    .95/(nrow+pad2) * (1-imscProp) * (1-padding)
                    ]');
                
                if ~isempty(freq)
                    plot(X, nanmean(newData, 2) );
                    set(handlesICA(fidx,2),'xlim',[X(1) X(end)])
                else
                    plot(X, nanmean(squeeze(data(idx, :, :)), 2));
                    set(handlesICA(fidx,2),'xlim',[myxlim(1) myxlim(2)])
                end
                set(handlesICA(fidx,2),'ylim',[-erpY  erpY])
                title(num2str(erpY), 'FontSize', 2 * erpY)
                %        ylim([-3 3])
                axis off
                %        colorbar
                %plotIdx = plotIdx + 1;
                %        subplot(nrows,ncols,[plotIdx+1+pFix plotIdx+1+ncols+pFix])
                handlesICA(fidx,3)=axes('position',...
                    [padding/2+(icol-1)/(ncol+pad2) + 1/(2*ncol)
                    padding/2+(nrow-irow)/(nrow+pad2)
                    .95/(ncol*2+pad2) * (1-padding*5)
                    .95/(nrow+pad2) * (1-padding)
                    ]');
                topoplot(EEG.icawinv(:, idx), EEG.chanlocs(EEG.icachansind));
                [~, kkkkk] = max( EEG.etc.ic_classification.ICLabel.classifications(idx, :) );
                tt = {[ num2str(idx), ' ', EEG.etc.ic_classification.ICLabel.classes{kkkkk}, ' ', num2str(100 * EEG.etc.ic_classification.ICLabel.classifications(idx, kkkkk), '%2.1f') ], ...
                        EEG.dipfit.model(idx).areadk };
                title(tt, 'Units', 'normalized', 'Position', [0.5, -0.1, 0])
                axis off
                
                %       if mod(idx,3) == 0, plotIdx=plotIdx+6; end
                
                %       [idx plotIdx plotIdx+ncols plotIdx+1+pFix plotIdx+1+ncols+pFix]
                
                
            end
        end
    end
elseif strcmpi(dataType,'data')
    dRange = mean(prctile(prctile(abs(data),95),95)); dRange=10;
    for idx=1:numel(EEG.chanlocs)
        %    try
        FC = 'k';
        if isempty(find(idx==EEG.icachansind))
            FC = 'r';
        end
        xp=sin(deg2rad(EEG.chanlocs(idx).theta))*EEG.chanlocs(idx).radius*positionScaleX+.5-figureSize/2;
        yp=cos(deg2rad(EEG.chanlocs(idx).theta))*EEG.chanlocs(idx).radius*positionScaleY+.5-figureSize/2+.05;
        
        if strcmpi({EEG.chanlocs(idx).labels},'PO5')
            xp=xp-.04;
        elseif strcmpi({EEG.chanlocs(idx).labels},'PO6')
            xp=xp+.04;
        elseif strcmpi({EEG.chanlocs(idx).labels},'PO7')
            xp=xp-.09;
        elseif strcmpi({EEG.chanlocs(idx).labels},'PO8')
            xp=xp+.09;
        elseif strcmpi({EEG.chanlocs(idx).labels},'P5')
            xp=xp-.02;
        elseif strcmpi({EEG.chanlocs(idx).labels},'P6')
            xp=xp+.02;
        elseif strcmpi({EEG.chanlocs(idx).labels},'P7')
            xp=xp-.06;
        elseif strcmpi({EEG.chanlocs(idx).labels},'P8')
            xp=xp+.06;
        elseif strcmpi({EEG.chanlocs(idx).labels},'O1')
            yp=yp-.03;
        elseif strcmpi({EEG.chanlocs(idx).labels},'Oz')
            yp=yp-.03;
        elseif strcmpi({EEG.chanlocs(idx).labels},'O2')
            yp=yp-.03;
        end
        
        handlesICA(idx)=axes('position',...
            [ xp
            yp
            figureSize
            figureSize
            ]');
        hold off
        
        if nargin == 4
            data=DATA;
            [wt,f]=cwt(reshape(data(idx,:,:),1,[]),EEG.srate);
            freqidx=dsearchn(f,freq');
            data=reshape(mean(wt(freqidx(2):freqidx(1),:)),EEG.srate*4-1,[]); % 4 = 4 seconds, -1 = 1 sample point
            data=abs(data);
            
            baseidx= dsearchn(EEG.times',[-400 0]');
            %            blc = mean(data(blcidx(1):blcidx(2),:));
            
            newData=10*log10( bsxfun(@rdivide, data(:,:), mean(data(baseidx(1):baseidx(2),:),1)));
            dRange = mean(prctile(abs(newData(idx,:)),95)); dRange=10;
            
            
            imagesc(X,(newData)')
            
            set(gca,'clim',[-dRange dRange])
            %title(num2str(dRange))
            axis off
            
            
        else
            imagesc(X,squeeze(data(idx,:,:))')
            set(gca,'clim',[-dRange dRange])
            %        set(gca,'XTickLabel',X)
            if FC == 'r', set(handlesICA(idx),'LineWidth',2,'Color',FC);
            else, axis off
            end
            title(EEG.chanlocs(idx).labels,'Color',FC)
            %axis off
        end
        
        
        %set(gca,'XLim',EEG.times(1:115:end))
        
        %    end
        
        handlesICA(idx)=axes('position',...
            [ xp
            yp-.02
            figureSize
            .02
            ]');
        hold off
        mdRange = max(abs(mean(squeeze(data(idx,:,:)),2)));
        mdRange=6;
        plot(mean(squeeze(data(idx,:,:)),2),'k-'), axis([xlim' -mdRange mdRange])
        
        axis off
        
        
        
        
    end
    axes('position',...
        [.95-figureSize/2
        .95-figureSize/2
        figureSize
        figureSize
        ]');
    hold off
    colorbar
    set(gca,'clim',[-dRange dRange])
    %        set(gca,'XTickLabel',X)
    axis off
    
end

end
