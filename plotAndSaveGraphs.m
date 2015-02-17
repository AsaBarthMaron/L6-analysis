function plotAnddoSaveGraphs(fileList,indices, channels, plotType,doSave)

if strcmp(plotType,'heatmaps')
    for k = indices
        path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
        load(path, 'outerSeq','innerSeq')
        [masterData,slaveData] = laserDelay_tuningCurve_firstAnalysis(path,'fast');
        
        for i = 1:length(masterData)
            diffData{i} = (masterData{i}-slaveData{i})./slaveData{i};
            %     diffData{i}(isnan(diffData{i})) = 0;
            %     diffData{i}(~isfinite(diffData{i})) = 0;
        end
        
        minRate = [min(min(masterData{1})),min(min(slaveData{1}))];
        maxRate = [max(max(masterData{1})),max(max(slaveData{1}))];
        minDiff = min(min(diffData{1}));
        maxDiff = max(max(diffData{1}));
        
        for i = 1:length(masterData)
            if minRate(1)>min(min(masterData{i}))
                minRate(1) = min(min(masterData{i}));
            end
            if minRate(2)>min(min(slaveData{i}))
                minRate(2) = min(min(slaveData{i}));
            end
            if minDiff(1)>min(min(diffData{i}))
                minDiff(1) = min(min(diffData{i}));
            end
            if maxRate(1)<max(max(masterData{i}))
                maxRate(1) = max(max(masterData{i}));
            end
            if maxRate(2)<max(max(slaveData{i}))
                maxRate(2) = max(max(slaveData{i}));
            end
            if maxDiff(1)<max(max(diffData{i}))
                maxDiff(1) = max(max(diffData{i}));
            end
        end
        
        figure
        j =1;
        rgbVals = 0:0.125:1;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.05 0.02], [0.05 0.02]);
        
        % Plot difference functions alongside heatmaps for each channel
        for i = [7,13,4,12,5,15,2,16,1,14,3,11,6]
            clims = [min(min([masterData{i},slaveData{i}])),max(max([masterData{i},slaveData{i}]))];
            subplot(13,3,j)
            imagesc(length(outerSeq.master.values),innerSeq.master.values,slaveData{i},clims)
            ylabel([num2str(i)])
            set(gca,'box','off')
            colormap(hot)
            axis xy
            set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
            j = j+1;
            
            subplot(13,3,j)
            imagesc(length(outerSeq.master.values),innerSeq.master.values,masterData{i},clims)
            set(gca,'box','off')
            colormap(hot)
            axis xy
            set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
            j = j+1;
            
            subplot(13,3,j)
            h = 1;
            for line = 2:2:length(diffData{i})
                plot(diffData{i}(line,:),'color', [rgbVals(h),0,1-rgbVals(h)],'linewidth',2)
                hold on
                h = h+1;
            end
            plot(1:8,zeros(8,1),'k--')
            set(gca,'box','off')
            % axis([outerSeq.master.values(1)./1000 outerSeq.master.values(8)./1000 min(min(diffData{i})) max(max(diffData{i}))])
            axis([1 8 -1 8])
            set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
            j = j+1;
            
        end
        subplot(13,3,1)
        title('Tone')
        subplot(13,3,2)
        title('Tone + Laser')
        subplot(13,3,3)
        title('Difference')
        subplot(13,3,37)
        xlabel('Frequency (KHz)')
        subplot(13,3,38)
        xlabel('Frequency (KHz)')
        subplot(13,3,39)
        xlabel('Frequency (KHz)')
        
        if doSave
            set(gcf,'position',[0,-1000,1200,1920]);
            export_fig(['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}(1:72) fileList{k}(73:end-4) '-Heatmaps.png'])
            close
        end
    end
elseif strcmp(plotType, 'heatmapDiff')
    Divisive = zeros(101,3);
    Divisive(1:50,3) = 1:-.02:.02;
    Divisive(52:end,1) = 0.02:.02:1;
    Divisive(52:end,2) = 0.02:.02:1;
    
    
    %        Divisive(1:50,3) = 1:-.01:0.51;
    %        Divisive(1:50,1) = 0.01:.01:0.5;
    %        Divisive(52:end,1) = .51:.01:1;
    %        Divisive(52:end,3) = 0.5:-.01:0.01;
    
    for k = indices
        path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
        load(path, 'outerSeq','innerSeq')
        [masterData,slaveData] = laserDelay_tuningCurve_firstAnalysis(path,'fast');
        
        
        figure
        j =1;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.05 0.02], [0.05 0.02]);
        
        % Plot difference functions alongside heatmaps for each channel
        for i = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6]
            subplot(8,2,j)
            imagesc(length(outerSeq.master.values),innerSeq.master.values,masterData{i}-slaveData{i})
            ylabel([num2str(i)])
            xlabel('Frequency (KHz)')
            set(gca,'box','off')
            colormap(Divisive)
            axis xy
            set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
            j = j+1;
            
            
        end
        subplot(8,2,1)
        title(fileList{k}(73:end))
        if doSave
            set(gcf,'position',[0,-1000,1200,1920]);
            %                export_fig(['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}(1:72) fileList{k}(73:end-4) '-HeatmapDiff.png'])
            export_fig(['C:\Users\polley_lab\Documents\MATLAB\Heatmaps\' fileList{k}(73:end-4) '-HeatmapDiff.png'])
            close
        end
    end
    
elseif strcmp(plotType, 'smoothPSTHs')
    scalingFactor = 1000/7;
    
    for k = indices
        path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
        load(path, 'outerSeq','innerSeq')
        [~,~,~,repPsthData,repPsthStdData] = laserDelay_tuningCurve_firstAnalysis(path,'smoothPSTH');
        channels = 4;
        
        % Pick a channel and average the PSTHs across repetitions
        repPSTHs = repPsthData{channels};
        repStdPSTHs = repPsthStdData{channels};
        % Average the PSTHs for a given channel and each delay of all
        % tones
        for i = 1:17
            for j = 1:8
                repNum(j) = size(repPSTHs{i,j},1);
            end
            %                maxRepNum = max(repNum);
            %                repAvgPSTHs{i} = zeros(maxRepNum, size(repPSTHs{i,1},2));
            %                for j = 1:8
            %                    repAvgPSTHs{i} = repAvgPSTHs{i}+padarray(repPSTHs{i,j},maxRepNum-repNum(j),'post');
            %                end
            minRepNum = min(repNum);
            repAvgPSTHs{i} = zeros(minRepNum, size(repPSTHs{i,1},2));
            for j = 1:8
                repAvgPSTHs{i} = repAvgPSTHs{i}+repPSTHs{i,j}(1:minRepNum,:);
            end
            repAvgPSTHs{i} = repAvgPSTHs{i}/j;
        end
        repAvgStdPSTHs = mean(mean(repStdPSTHs));
        
        figure
        spacer = 3*repAvgStdPSTHs*scalingFactor;
        hold on
        for i =1:17
            alpha = .8;
            plot(squeeze(repAvgPSTHs{i}*scalingFactor)'+(i*spacer),'color',[0 0 0]+alpha)
            plot(mean(repAvgPSTHs{i}*scalingFactor)+(i*spacer),'linewidth',2,'color',[0 0 0])
            h = patch([(150+(i*50)) (150+(i*50)) (550+(i*50)) (550+(i*50))],[((-spacer/6)+(i*spacer)) ((spacer/1.4)+(i*spacer)) ((spacer/1.4)+(i*spacer)) ((-spacer/6)+(i*spacer))],'b');
            set(h, 'EdgeColor','none')
            set(h, 'FaceAlpha',.3)
        end
        y = get(gca,'YLim')
        plot(1000*ones(y(2)+1,1),(y(1):y(2)),'--','linewidth',2)
        plot(3000*ones(y(2)+1,1),(y(1):y(2)),'--r','linewidth',2)
        title(['Channel ' num2str(channels)])
        hold off
        if doSave
            set(gcf,'position',[0,-1000,1200,1920]);
            %                export_fig(['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}(1:72) fileList{k}(73:end-4) '-HeatmapDiff.png'])
            export_fig(['C:\Users\polley_lab\Documents\MATLAB\PSTHs\' fileList{k}(73:end-4) '-PSTHs.png'])
            close
        end
    end
    
elseif strcmp(plotType, 'PSTHs')
    
    for k = indices
        path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
        load(path, 'outerSeq','innerSeq')
        [~,~,~,repPsthData,repPsthStdData] = laserDelay_tuningCurve_firstAnalysis(path,'PSTH');
        for ch = channels';
            
            % Pick a channel and average the PSTHs across repetitions
            repPSTHs = repPsthData{ch};
            repStdPSTHs = repPsthStdData{ch};
            % Average the PSTHs for a given channel and each delay of all
            % tones
            for i = 1:17
                for j = 1:8
                    repNum(j) = size(repPSTHs{i,j},1);
                end
                %                maxRepNum = max(repNum);
                %                repAvgPSTHs{i} = zeros(maxRepNum, size(repPSTHs{i,1},2));
                %                for j = 1:8
                %                    repAvgPSTHs{i} = repAvgPSTHs{i}+padarray(repPSTHs{i,j},maxRepNum-repNum(j),'post');
                %                end
                % Exclude reps if they do not have a repeition for each
                % trial combination
                minRepNum = min(repNum);
                repAvgPSTHs{i} = zeros(minRepNum, size(repPSTHs{i,1},2));
                for j = 1:8
                    repAvgPSTHs{i} = repAvgPSTHs{i}+repPSTHs{i,j}(1:minRepNum,:);
                end
                repAvgPSTHs{i} = repAvgPSTHs{i}/j;
            end
            repAvgStdPSTHs = mean(mean(repStdPSTHs));
            
            figure
            spacer = 2*repAvgStdPSTHs;
            hold on
            for i =1:17
                meanRepAvgPSTHs = mean(repAvgPSTHs{i});
                % Find spontaneous FR std (at 2000ms - bin 400)
                spntStd = std(meanRepAvgPSTHs(400:499));
                spntMean = mean(meanRepAvgPSTHs(400:499));
                tmp = zeros(1,3500);
                tmp(1:5:3500) = meanRepAvgPSTHs;
                tmp(2:5:3500) = meanRepAvgPSTHs;
                tmp(3:5:3500) = meanRepAvgPSTHs;
                tmp(4:5:3500) = meanRepAvgPSTHs;
                tmp(5:5:3500) = meanRepAvgPSTHs;
                meanTrace = tmp; clear tmp;
                %          tmp = zeros(minRepNum,1,3500);
                %           tmp(:,1:5:3500) = repAvgPSTHs{i}(1:minRepNum,:);
                %            tmp(:,2:5:3500) = repAvgPSTHs{i}(1:minRepNum,:);
                %             tmp(:,3:5:3500) = repAvgPSTHs{i}(1:minRepNum,:);
                %              tmp(:,4:5:3500) = repAvgPSTHs{i}(1:minRepNum,:);
                %               tmp(:,5:5:3500) = repAvgPSTHs{i}(1:minRepNum,:);
                %                traces = tmp; clear tmp;
                alpha = .8;
                %                plot(1:5:3500,meanRepAvgPSTHs+(i*spacer),'k','linewidth',2)
                %               plot(squeeze(traces)'+(i*spacer),'color',[0 0 0]+alpha)
                plot(meanTrace+(i*spacer),'k','linewidth',2)
                plot(1:3500,3*spntStd*ones(3500,1)+(i*spacer)+spntMean,'k--')
                h = patch([(150+(i*50)) (150+(i*50)) (550+(i*50)) (550+(i*50))],[((-spacer/6)+(i*spacer)) ((spacer/1.4)+(i*spacer)) ((spacer/1.4)+(i*spacer)) ((-spacer/6)+(i*spacer))],'b');
                set(h, 'EdgeColor','none')
                set(h, 'FaceAlpha',.3)
            end
            y = get(gca,'YLim')
            plot(1000*ones(y(2)+1,1),(y(1):y(2)),'--','linewidth',2)
            plot(3000*ones(y(2)+1,1),(y(1):y(2)),'--r','linewidth',2)
            title(['Experiment: ' fileList{k}(73:end-4) ' - Channel ' num2str(ch) ])
            hold off
            if doSave
                if ~isdir(['C:\Users\polley_lab\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\' fileList{k}(73:98) '\keep\'])
                    mkdir(['C:\Users\polley_lab\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\' fileList{k}(73:98) '\keep\']);
                end
                set(gcf,'position',[0,-1000,1200,1920]);
                export_fig(['C:\Users\polley_lab\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Images\' fileList{k}(73:98) '\keep\' fileList{k}(73:end-4) '-channel_' num2str(ch) '.png'])
                close
            end
        end
    end
end