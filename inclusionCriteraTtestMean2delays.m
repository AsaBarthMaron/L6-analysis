function [summaryResponses,pValsChansExp] = inclusionCriteraTtest(fileList,indices)

% for k = indices
%     path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
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
    path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
    load(path, 'outerSeq','innerSeq')
    [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(path,'PSTH',1);
    channels = 1:16;
    
    for i = channels;
        meanPSTHs = squeeze(mean(psthData{i},2));
        spntFR = meanPSTHs(:,2970:2999);
        toneFR = meanPSTHs(:,3016:3045);
        
        j=1;
        for h = 1:2:16
            [H,P] = ttest2(mean([spntFR(h,:);spntFR(h+1,:)]),mean([toneFR(h,:);toneFR(h+1,:)]),.05,'left');
            toneResponsive(j) = H;
            pVals(j) = P;
            j = j+1;
            clear P H;
        end
        [H,P] = ttest2(mean([spntFR(16,:);spntFR(16+1,:)]),mean([toneFR(16,:);toneFR(16+1,:)]),.05,'left');
        toneResponsive(9) = H;
        pVals(9) = P;   
        clear P H;
        
        toneResponsiveChans(i,:) = toneResponsive;
        pValsChans(i,:) = pVals;
        clear toneResponsive pVals;
    end
    toneResponsesChansExp{k} = toneResponsiveChans;
    toneResponsesExp = sum(toneResponsiveChans')';
    summaryResponses{k} = find(toneResponsesExp>=9);
    pValsChansExp{k} = pValsChans;
    clear toneResponsiveChans pValsChans;
end

