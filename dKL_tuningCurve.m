L23masterData = [];
L4masterData = [];
L5masterData = [];
L6masterData = [];
L23slaveData = [];
L4slaveData = [];
L5slaveData = [];
L6slaveData = [];
L23psthData = [];
L4psthData = [];
L5psthData = [];
L6psthData = [];
% 
%  % Rectified  
% for k = indices
%     
%     %     for i = toneResponses{k}
%     %     	gainQuant{k,i} = (mean(slaveDataExp)-mean(masterDataExp{}))/(mean(slaveDataExp)+mean(masterDataExp{}))
%     %     end
%     % First reorganize all the data from exp>channel format to layer>exp
%     % format
%     for ch = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6]
%         
%         if (max(channelAssignment(k,1:3) == ch)) && (max(toneResponses{k} == ch))
%             % Layer 2/3
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster<0) = 0; rectMaster = rectMaster +.5;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave<0) = 0; rectSlave = rectSlave +.5;
%             L23masterData = [L23masterData {rectMaster}];
%             L23slaveData = [L23slaveData {rectSlave}];
%             L23psthData = [L23psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster<0) = 0; rectMaster = rectMaster +.5;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave<0) = 0; rectSlave = rectSlave +.5;
%             L4masterData = [L4masterData {rectMaster}];
%             L4slaveData = [L4slaveData {rectSlave}];
%             L4psthData = [L4psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster<0) = 0; rectMaster = rectMaster +.5;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave<0) = 0; rectSlave = rectSlave +.5;
%             L5masterData = [L5masterData {rectMaster}];
%             L5slaveData = [L5slaveData {rectSlave}];
%             L5psthData = [L5psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster<0) = 0; rectMaster = rectMaster +.5;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave<0) = 0; rectSlave = rectSlave +.5;
%             L6masterData = [L6masterData {rectMaster}];
%             L6slaveData = [L6slaveData {rectSlave}];
%             L6psthData = [L6psthData psthData(k,ch)];
%         end
%     end
% end

% No manipulation
for k = indices
    
    %     for i = toneResponses{k}
    %     	gainQuant{k,i} = (mean(slaveDataExp)-mean(masterDataExp{}))/(mean(slaveDataExp)+mean(masterDataExp{}))
    %     end
    % First reorganize all the data from exp>channel format to layer>exp
    % format
    for ch = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6]
        
        if (max(channelAssignment(k,1:3) == ch)) && (max(toneResponses{k} == ch))
            % Layer 2/3
            L23masterData = [L23masterData {masterDataExp{k,ch}+.5}];
            L23slaveData = [L23slaveData {slaveDataExp{k,ch}+.5}];
            L23psthData = [L23psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
            L4masterData = [L4masterData {masterDataExp{k,ch}+.5}];
            L4slaveData = [L4slaveData {slaveDataExp{k,ch}+.5}];
            L4psthData = [L4psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
            L5masterData = [L5masterData {masterDataExp{k,ch}+.5}];
            L5slaveData = [L5slaveData {slaveDataExp{k,ch}+.5}];
            L5psthData = [L5psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
            L6masterData = [L6masterData {masterDataExp{k,ch}+.5}];
            L6slaveData = [L6slaveData {slaveDataExp{k,ch}+.5}];
            L6psthData = [L6psthData psthData(k,ch)];
        end
    end
 end

masterData = [{L23masterData} {L4masterData} {L5masterData} {L6masterData}];
slaveData = [{L23slaveData} {L4slaveData} {L5slaveData} {L6slaveData}];
numTones = size(masterData{1}{1},2);

for layerIndex = 1:length(masterData)
    for muIndex = 1:length(masterData{layerIndex})
        numTones = size(masterData{layerIndex}{muIndex},2);
        % Calculate T(d,f) and TLdelay(d,f), f = frequency d = delay
        masterData{layerIndex}{muIndex} = masterData{layerIndex}{muIndex}-repmat(min(masterData{layerIndex}{muIndex},[],2),1,numTones)+0.5;
        slaveData{layerIndex}{muIndex} = slaveData{layerIndex}{muIndex}-repmat(min(slaveData{layerIndex}{muIndex},[],2),1,numTones)+0.5;
        TL = masterData{layerIndex}{muIndex}./repmat(sum(masterData{layerIndex}{muIndex},2),1,numTones);
        T = slaveData{layerIndex}{muIndex}./repmat(sum(slaveData{layerIndex}{muIndex},2),1,numTones);
        dKL{layerIndex}(:,muIndex) = diag(T * log2(T./TL)');
        
        % control measurement taken by calculating dKL between a control tuning curve at
        % a given delay to control tuning curves at all other delays
        for delayIndex = 1:size(masterData{layerIndex}{muIndex},1)
            dKLnoiseTmp = diag(T * log2(T./repmat(T(delayIndex,:),17,1))');
            dKLnoiseTmp(delayIndex) = nan;
            dKLnoise{layerIndex}(delayIndex,muIndex) = nanmean(dKLnoiseTmp);
        end
%         % control measurement taken by calculating all pairwise
%         % combinations of control tuning curves for a given MU
%         delayCombinations = nchoosek(1:size(masterData{layerIndex}{muIndex},1),2);
%         for combIndex = 1:size(delayCombinations,1)
%             dKLnoiseTmp(2*combIndex) = T(delayCombinations(combIndex,1),:) * ...
%                 log2(T(delayCombinations(combIndex,1),:)./T(delayCombinations(combIndex,2),:))';
%             dKLnoiseTmp(2*combIndex-1) = T(delayCombinations(combIndex,2),:) * ...
%                 log2(T(delayCombinations(combIndex,2),:)./T(delayCombinations(combIndex,1),:))';
%         end
%         dKLnoise{layerIndex}(:,muIndex) = ones(17,1).* mean(dKLnoiseTmp);
    end
end

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
layerName = [{'Layer 2/3'} {'Layer 4'} {'Layer 5'} {'Layer 6'} ];
for layerIndex = 1:length(masterData)
    meanDKL = fliplr(mean(dKL{layerIndex},2)');
    semDKL = fliplr(std(dKL{layerIndex},0,2)'./sqrt(length(masterData{layerIndex})));
    meanDKLnoise = fliplr(mean(dKLnoise{layerIndex},2)');
    semDKLnoise = fliplr(std(dKLnoise{layerIndex},0,2)'./sqrt(length(masterData{layerIndex})));
    subplot(2,2,layerIndex)
    area([0:25:400,400:-25:0],[meanDKLnoise-semDKLnoise,fliplr(meanDKLnoise+semDKLnoise)],'facecolor',[0 .7 .9],'linestyle','none');
    hold on
    area([0:25:400,400:-25:0],[meanDKL-semDKL,fliplr(meanDKL+semDKL)],'facecolor',[.9 0.1 0.1],'linestyle','none');
    plot(0:25:400,meanDKLnoise,'color',[0 .2 .3],'linewidth',3)
    plot(0:25:400,meanDKL,'color',[0.3 0 0],'linewidth',3)
    %axis([0 400 0 1.5])
    title(layerName{layerIndex})
    xlabel('Tone onset (relative to laser onset)')
    ylabel('Gain')
    set(gca,'box','off')
end
    legend('\DeltaH control','\DeltaH laser')
    legend('boxoff')

