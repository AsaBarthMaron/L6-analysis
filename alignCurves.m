for layerIndex = 1:4
    numMUs = length(masterData{layerIndex});
    numTones = size(masterData{layerIndex}{1},2);
    numSlots = (2*numTones)-1;
    numDelays = size(masterData{layerIndex}{1},1);
    
    alignedMaster{layerIndex} = nan(numDelays,numSlots,numMUs);
    alignedSlave{layerIndex} = nan(numDelays,numSlots,numMUs);
    
    for muIndex = 1:numMUs
        for delayIndex = 1: numDelays
           tuningCurve = masterData{layerIndex}{muIndex}(delayIndex,:);
           bfIndex = find(tuningCurve == max(tuningCurve));
           centeredRange = (numTones-(bfIndex-1)):(numSlots-(bfIndex-1));
           alignedMaster{layerIndex}(delayIndex,centeredRange,muIndex) = tuningCurve;
           
           tuningCurve = slaveData{layerIndex}{muIndex}(delayIndex,:);
           bfIndex = find(tuningCurve == max(tuningCurve));
           centeredRange = (numTones-(bfIndex-1)):(numSlots-(bfIndex-1));
           alignedSlave{layerIndex}(delayIndex,centeredRange,muIndex) = tuningCurve;
        end
    end
    meanAlMaster{layerIndex} = nanmean(alignedMaster{layerIndex},3);
    meanAlSlave{layerIndex} = nanmean(alignedSlave{layerIndex},3);
    
    numCurves = sum(~isnan(alignedMaster{layerIndex}),3);
    numCurves(numCurves == 0) = nan;
    semAlMaster{layerIndex} = nanstd(alignedMaster{layerIndex},0,3)./sqrt(numCurves);
    semAlMaster{layerIndex}(isnan(semAlMaster{layerIndex})) = 0;
    
    numCurves = sum(~isnan(alignedSlave{layerIndex}),3);
    numCurves(numCurves == 0) = nan;
    semAlSlave{layerIndex} = nanstd(alignedSlave{layerIndex},0,3)./sqrt(numCurves);
    semAlSlave{layerIndex}(isnan(semAlSlave{layerIndex})) = 0;
end

for layerIndex = 1:4
    figure
    subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.03 0.03], [0.05 0.02]);
    relLaserIndex = numDelays;
    
    for delayIndex = 1:numDelays
        masterSlice = meanAlMaster{layerIndex}(relLaserIndex,:);
        slaveSlice = meanAlSlave{layerIndex}(relLaserIndex,:);
        masterSlice(isnan(masterSlice)) = 0;
        slaveSlice(isnan(slaveSlice)) = 0;
        
        subplot(9,2,delayIndex)
        %plot(squeeze(alignedSlave{layerIndex}(delayIndex,:,:)),'color',[0 .9 0])
        hold on
        %plot(squeeze(alignedMaster{layerIndex}(delayIndex,:,:)),'color',[0 0 .9])
        
        area([1:15,15:-1:1],[slaveSlice-semAlSlave{layerIndex}(relLaserIndex,:),...
            fliplr(slaveSlice+semAlSlave{layerIndex}(relLaserIndex,:))],...
            'facecolor',[ 0 .9 0],'linestyle','none');
        area([1:15,15:-1:1],[masterSlice-semAlMaster{layerIndex}(relLaserIndex,:),...
            fliplr(masterSlice+semAlMaster{layerIndex}(relLaserIndex,:))],...
            'facecolor',[ 0 0 .9],'linestyle','none');
        
        plot(meanAlSlave{layerIndex}(relLaserIndex,:),'color',[0 .3 0],'linewidth',3)
        plot(meanAlMaster{layerIndex}(relLaserIndex,:),'color',[0 0 .3],'linewidth',3)
        plot(1:numSlots,zeros(numSlots,1),'k--')
        ylabel([num2str((delayIndex*25)-25) ' ms'])
        
        subplot(9,2,1)
            if layerIndex == 1
                title('Layer 23');
            elseif layerIndex == 2
                title('Layer 4');
            elseif layerIndex == 3
                title('Layer 5');
            elseif layerIndex == 4
                title('Layer 6');
            end
        relLaserIndex = relLaserIndex-1;
    end
end

