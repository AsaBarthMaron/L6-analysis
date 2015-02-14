function summaryResponses = inclusionCriteria1(fileList,indices)
% For each experiment passed to the function
for k = indices
    path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
    load(path, 'outerSeq','innerSeq')
    [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(path,'PSTH');
    channels = 1:16;

    % For each of the 16 channels
    for i = channels;
        meanPSTHs = squeeze(mean(psthData{i},2));

        spntFR = mean(meanPSTHs(:,400:499),2);
        spntStd = std(meanPSTHs(:,400:499)');
        toneResponsive = [];
        % For each delay on the channel
        for j = 1:length(innerSeq.master.values)
            points = 0;
            firstBin = [];
            startTracking = 0;
            win = 0;
            for johnnieWalker = 600:609
                if meanPSTHs(j,johnnieWalker)>((3*spntStd(j))+spntFR(j))
                    points = points + 1;
                    if points == 1
                        firstBin = johnnieWalker;
                        startTracking = 1;
                    end     
                end
                if (meanPSTHs(j,johnnieWalker)<=((3*spntStd(j))+spntFR(j))) && startTracking
                    points = 0;
                    firstBin = [];
                    startTracking = 0;
                end
                if points >= 2
                    win = 1;
                    break
                end
            end
            if win && (mean(meanPSTHs(j,firstBin:firstBin+4))> ((3*spntStd(j))+spntFR(j)))
                toneResponsive(j) = 1;
            else
                toneResponsive(j) = 0;
            end
        end
        toneResponsiveChans(i,:) = toneResponsive;
    end
    toneResponsesChansExp{k} = toneResponsiveChans;
    toneResponsesExp = sum(toneResponsiveChans')';
    summaryResponses{k} = find(toneResponsesExp==17);
    clear toneResponsiveChans;
end
            