function [summaryResponses,pValsChansExp] = inclusionCriteraTtest(fileList,indices)

% for k = indices
%     path = ['C:\Users\asa\Documents\MATLAB\' fileList{k}];
%     load(path, 'outerSeq','innerSeq')
%     [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(path,'PSTH',1);
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
    path = ['C:\Users\asa\Documents\MATLAB\' fileList{k}];
    load(path, 'outerSeq','innerSeq')
    [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(path,'PSTH',1);
    channels = 1:16;
    
    for i = channels;
        meanPSTHs = squeeze(mean(psthData{i},2));
        spntFR = meanPSTHs(:,2970:2999);
        toneFR = meanPSTHs(:,3016:3045);
        
        for j = 1:length(innerSeq.master.values)
            [H,P] = ttest2(spntFR(j,:),toneFR(j,:),'Alpha',.01,'tail','left');
            toneResponsive(j) = H;
            pVals(j) = P;
            clear P H;
        end
        toneResponsiveChans(i,:) = toneResponsive;
        pValsChans(i,:) = pVals;
        clear toneResponsive pVals;
    end
    toneResponsesChansExp{k} = toneResponsiveChans;
    toneResponsesExp = sum(toneResponsiveChans')';
    summaryResponses{k} = find(toneResponsesExp>=17);
    pValsChansExp{k} = pValsChans;
    clear toneResponsiveChans pValsChans;
end

