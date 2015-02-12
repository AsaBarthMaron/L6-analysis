        scalingFactor = 1000/7;
        
        for k = indices
            path = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
            load(path, 'outerSeq','innerSeq')
            [~,~,psthData] = laserDelay_tuningCurve_firstAnalysis(path,'deep');
            channels = 13;
            
            for i = channels;
            	meanPSTHs = squeeze(mean(psthData{i},2))*scalingFactor;
                
                spntFR = mean(meanPSTHs(:,2000:2020),2);
                for j = 1:length(innerSeq.master.values)
                    spntFR(j) = mean(meanPSTHs(i,2000:2020));
                end
            end
        end
            
                