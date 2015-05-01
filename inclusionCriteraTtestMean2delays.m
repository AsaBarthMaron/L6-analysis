function [summaryResponses,pValsChansExp] = inclusionCriteraTtestMean2delays(fileList,indices)

% for k = indices
%     filePath = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
%     load(filePath, 'outerSeq','innerSeq')
%     [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(filePath,'PSTH',1);
%     channels = 1:16;
%     
%     for i = channels;
%         meanPSTHs = squeeze(mean(psthData{i},2));
%         spntFR = mean(meanPSTHs(:,2970:2999));
%         toneFR = mean(meanPSTHs(:,3016:3045));
%         
%         [H,P] = ttest2(spntFR,toneFR,'Alpha',.01,'tail','left');
%         toneResponsiveChans(i) = H;
%         pValsChans(i) = P;
%         clear P H spntFR toneFR;
% 
%     end
%     toneResponsesChansExp{k} = toneResponsiveChans;
%     pValsChansExp{k} = pValsChans;
%     clear toneResponsiveChans pValsChans;
% end
%  

for k = indices
    filePath{k} = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
end
psthDataExps = {};
parfor k = indices
    [~,~,psthDataExps{k}] = laserDelay_tuningCurve_firstAnalysis(filePath{k},'PSTH',1);
end
for k = indices
    psthData = psthDataExps{k};
    load(filePath{k}, 'outerSeq','innerSeq','stimChans')
    channels = 1:16;
    
    for i = channels;
        meanPSTHs = squeeze(mean(psthData{i},2));
        spntFR = meanPSTHs(:,stimChans{1,4}.Gate.Delay-30:stimChans{1,4}.Gate.Delay-1);
        toneFR = meanPSTHs(:,stimChans{1,4}.Gate.Delay+16:stimChans{1,4}.Gate.Delay+45);
        
        j=1;
        for h = 1:2:(size(meanPSTHs,1)-1)
            [H,P] = ttest2(mean([spntFR(h,:);spntFR(h+1,:)]),mean([toneFR(h,:);toneFR(h+1,:)]),.01,'left');
            toneResponsive(j) = H;
            pVals(j) = P;
            j = j+1;
            clear P H;
        end
        if mod(size(meanPSTHs,1),2)
            [H,P] = ttest2(mean([spntFR(size(meanPSTHs,1)-1,:);spntFR(size(meanPSTHs,1),:)]),...
                mean([toneFR(size(meanPSTHs,1)-1,:);toneFR(size(meanPSTHs,1),:)]),.01,'left');
            toneResponsive(size(meanPSTHs,1)) = H;
            pVals(size(meanPSTHs,1)) = P;   
            clear P H;
        end
        toneResponsiveChans(i,:) = toneResponsive;
        pValsChans(i,:) = pVals;
        clear toneResponsive pVals;
    end
    toneResponsesChansExp{k} = toneResponsiveChans;
    toneResponsesExp = sum(toneResponsiveChans')';
    summaryResponses{k} = find(toneResponsesExp>=(ceil(size(meanPSTHs,1)/2)));
    pValsChansExp{k} = pValsChans;
    clear toneResponsiveChans pValsChans;
end

