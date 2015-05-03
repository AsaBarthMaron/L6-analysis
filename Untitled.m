[masterData, slaveData] = laserDelay_tuningCurve_firstAnalysis('C:\Users\asa\Documents\MATLAB\cortex_genericlaser+tuningmoddelay\Experimental\Ntsr1-cre#037\asa030615\asa030615-1-21-tuningcurve.mat','fast',5);

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
        
figure
imagesc(masterData{15})
colormap(hot)
set(gca,'YTick',[1:18]);
set(gca,'YTickLabel',{'50(2.5ms)','50(5ms)','50(10ms)','50(25ms)','50(50ms)',	'50(100ms)'	,'50(200ms)','50(300ms)','50(400ms)',...
    '150(2.5ms)','150(5ms)','150(10ms)','150(25ms)','150(50ms)','150(100ms)','150(200ms)','150(300ms)','150(400ms)'});
