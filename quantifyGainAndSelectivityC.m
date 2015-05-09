%function[L23masterData, L4masterData, L5masterData, L6masterData] = quantifyGain(indices,toneResponses,startingBins,channelAssignment) 
%for k = indices
%    filePath{k} = ['C:\Users\polley_lab\Documents\MATLAB\' fileList{k}];
%end
%psthData = {};
%parfor k = indices
%    [~,~,psthData(k,:)] = laserDelay_tuningCurve_firstAnalysis(filePath{k},'smoothPSTH',2);
%end    
%[masterDataExp,slaveDataExp] = laserDelay_tuningCurve_secondAnalysisTtest(indices,toneResponses,0,fileList)
%close all

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
            L23masterData = [L23masterData masterDataExp(k,ch)];
            L23slaveData = [L23slaveData slaveDataExp(k,ch)];
            L23psthData = [L23psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
            L4masterData = [L4masterData masterDataExp(k,ch)];
            L4slaveData = [L4slaveData slaveDataExp(k,ch)];
            L4psthData = [L4psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
            L5masterData = [L5masterData masterDataExp(k,ch)];
            L5slaveData = [L5slaveData slaveDataExp(k,ch)];
            L5psthData = [L5psthData psthData(k,ch)];
        elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
            L6masterData = [L6masterData masterDataExp(k,ch)];
            L6slaveData = [L6slaveData slaveDataExp(k,ch)];
            L6psthData = [L6psthData psthData(k,ch)];
        end
    end
end

%  %Zero shifted down.
% for k = indices
%     
% %         for i = toneResponses{k}
% %         	gainQuant{k,i} = (mean(slaveDataExp)-mean(masterDataExp{}))/(mean(slaveDataExp)+mean(masterDataExp{}))
% %         end
% %     First reorganize all the data from exp>channel format to layer>exp
% %     format
%     for ch = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6]
%         
%         if (max(channelAssignment(k,1:3) == ch)) && (max(toneResponses{k} == ch))
%             % Layer 2/3
%             zeroVal = min(min([masterDataExp{k,ch} slaveDataExp{k,ch}]));
%             L23masterData = [L23masterData {masterDataExp{k,ch}-zeroVal+.5}];
%             L23slaveData = [L23slaveData {slaveDataExp{k,ch}-zeroVal+.5}];
%             L23psthData = [L23psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
%             zeroVal = min(min([masterDataExp{k,ch} slaveDataExp{k,ch}]));
%             L4masterData = [L4masterData {masterDataExp{k,ch}-zeroVal+.5}];
%             L4slaveData = [L4slaveData {slaveDataExp{k,ch}-zeroVal+.5}];
%             L4psthData = [L4psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
%             zeroVal = min(min([masterDataExp{k,ch} slaveDataExp{k,ch}]));
%             L5masterData = [L5masterData {masterDataExp{k,ch}-zeroVal+.5}];
%             L5slaveData = [L5slaveData {slaveDataExp{k,ch}-zeroVal+.5}];
%             L5psthData = [L5psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
%             zeroVal = min(min([masterDataExp{k,ch} slaveDataExp{k,ch}]));
%             L6masterData = [L6masterData {masterDataExp{k,ch}-zeroVal+.5}];
%             L6slaveData = [L6slaveData {slaveDataExp{k,ch}-zeroVal+.5}];
%             L6psthData = [L6psthData psthData(k,ch)];
%         end
%     end
% end

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

% % Rectified below zero
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
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster>0) = 0; rectMaster = abs(rectMaster)+.1;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave>0) = 0; rectSlave = abs(rectSlave)+.1;
%             L23masterData = [L23masterData {rectMaster}];
%             L23slaveData = [L23slaveData {rectSlave}];
%             L23psthData = [L23psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster>0) = 0; rectMaster = abs(rectMaster)+.1;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave>0) = 0; rectSlave = abs(rectSlave)+.1;
%             L4masterData = [L4masterData {rectMaster}];
%             L4slaveData = [L4slaveData {rectSlave}];
%             L4psthData = [L4psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster>0) = 0; rectMaster = abs(rectMaster)+.1;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave>0) = 0; rectSlave = abs(rectSlave)+.1;
%             L5masterData = [L5masterData {rectMaster}];
%             L5slaveData = [L5slaveData {rectSlave}];
%             L5psthData = [L5psthData psthData(k,ch)];
%         elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
%             rectMaster = masterDataExp{k,ch}; rectMaster(rectMaster>0) = 0; rectMaster = abs(rectMaster)+.1;
%             rectSlave = slaveDataExp{k,ch}; rectSlave(rectSlave>0) = 0; rectSlave = abs(rectSlave)+.1;
%             L6masterData = [L6masterData {rectMaster}];
%             L6slaveData = [L6slaveData {rectSlave}];
%             L6psthData = [L6psthData psthData(k,ch)];
%         end
%     end
% end
%% Gain quantification

for i =1:length(L23masterData)
    L23gainData(:,i) = ((mean(L23masterData{i},2)- mean(L23slaveData{i},2))...
                      ./mean(L23slaveData{i},2))*100;
    meanToneL23masterData(:,i) = mean(L23masterData{i},2);
    meanToneL23slaveData(:,i) = mean(L23slaveData{i},2);
end

for i =1:length(L4masterData)
    L4gainData(:,i) = ((mean(L4masterData{i},2)- mean(L4slaveData{i},2))...
                      ./mean(L4slaveData{i},2))*100;
    meanToneL4masterData(:,i) = mean(L4masterData{i},2);
    meanToneL4slaveData(:,i) = mean(L4slaveData{i},2);
end

for i =1:length(L5masterData)
    L5gainData(:,i) = ((mean(L5masterData{i},2)- mean(L5slaveData{i},2))...
                      ./mean(L5slaveData{i},2))*100;
    meanToneL5masterData(:,i) = mean(L5masterData{i},2);
    meanToneL5slaveData(:,i) = mean(L5slaveData{i},2);
end

for i =1:length(L6masterData)
    L6gainData(:,i) = ((mean(L6masterData{i},2)- mean(L6slaveData{i},2))...
                      ./mean(L6slaveData{i},2))*100;
    meanToneL6masterData(:,i) = mean(L6masterData{i},2);
    meanToneL6slaveData(:,i) = mean(L6slaveData{i},2);
end


%% Calculate Laser PSTH
L23psth = [];
L4psth = [];
L5psth = [];
L6psth = [];

for i =1:length(L23masterData)
    L23psth(i,:) = squeeze(mean(L23psthData{i}(1,:,:),2));
end
L23psth = (L23psth-mean(mean(L23psth(:,1000:1190))))./length(L23masterData);

for i =1:length(L4masterData)
    L4psth(i,:) = squeeze(mean(L4psthData{i}(1,:,:),2));
end
L4psth = (L4psth-mean(mean(L4psth(:,1000:1190))))./length(L4masterData);

for i =1:length(L5masterData)
    L5psth(i,:) = squeeze(mean(L5psthData{i}(1,:,:),2));
end
L5psth = (L5psth-mean(mean(L5psth(:,1000:1190))))./length(L5masterData);

for i =1:length(L6masterData)
   L6psth(i,:) = squeeze(mean(L6psthData{i}(1,:,:),2));
end
L6psth = (L6psth-mean(mean(L6psth(:,1000:1190))))./length(L6masterData);

maxPsth = max(max([mean(L23psth);mean(L4psth);mean(L5psth);mean(L6psth)]));
L23psth = L23psth./maxPsth;
L4psth = L4psth./maxPsth;
L5psth = L5psth./maxPsth;
L6psth =  L6psth./maxPsth;

%% Calculate ALDR/RLDR
for i = 1:length(L23masterData)
%     tmp = find(nanmean(L23masterData{i},1) == max(nanmean(L23masterData{i},1)));
%     BFL23masterData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L23masterData{i},1)
        normL23masterData{i}(j,:) = (L23masterData{i}(j,:)./max(L23masterData{i}(j,:)));
        bfIndex = find(normL23masterData{i}(j,:) == 1);
        L23masterALDR(j,i) = normL23masterData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L23masterALDR(j,i) = (L23masterALDR(j,i)+normL23masterData{i}(j,bfIndex(1)+1))./2;
            L23masterALDR(j,i) = L23masterALDR(j,i)./nanmean(normL23masterData{i}(j,[bfIndex(1)+2:size(normL23masterData{i},2)]));
        elseif (bfIndex(1)+1)> size(L23masterData{i},2)
            L23masterALDR(j,i) = (L23masterALDR(j,i)+normL23masterData{i}(j,bfIndex(1)-1))./2;
            L23masterALDR(j,i) = L23masterALDR(j,i)./nanmean(normL23masterData{i}(j,[1:bfIndex(1)-2]));
        elseif normL23masterData{i}(j,bfIndex(1)-1)>= normL23masterData{i}(j,bfIndex(1)+1)
            L23masterALDR(j,i) = (L23masterALDR(j,i)+normL23masterData{i}(j,bfIndex(1)-1))./2;
            L23masterALDR(j,i) = L23masterALDR(j,i)./nanmean(normL23masterData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL23masterData{i},2)]));
        elseif normL23masterData{i}(j,bfIndex(1)-1)< normL23masterData{i}(j,bfIndex(1)+1)
            L23masterALDR(j,i) = (L23masterALDR(j,i)+normL23masterData{i}(j,bfIndex(1)+1))./2;
            L23masterALDR(j,i) = L23masterALDR(j,i)./nanmean(normL23masterData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL23masterData{i},2)]));
        end
    end
    
end
for i = 1:length(L23slaveData)
%     tmp = find(nanmean(L23slaveData{i},1) == max(nanmean(L23slaveData{i},1)));
%     BFL23slaveData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L23slaveData{i},1)
        normL23slaveData{i}(j,:) = (L23slaveData{i}(j,:)./max(L23slaveData{i}(j,:)));
        bfIndex = find(normL23slaveData{i}(j,:) == 1);
        L23slaveALDR(j,i) = normL23slaveData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L23slaveALDR(j,i) = (L23slaveALDR(j,i)+normL23slaveData{i}(j,bfIndex(1)+1))./2;
            L23slaveALDR(j,i) = L23slaveALDR(j,i)./nanmean(normL23slaveData{i}(j,[bfIndex(1)+2:size(normL23slaveData{i},2)]));
        elseif (bfIndex(1)+1)> size(L23slaveData{i},2)
            L23slaveALDR(j,i) = (L23slaveALDR(j,i)+normL23slaveData{i}(j,bfIndex(1)-1))./2;
            L23slaveALDR(j,i) = L23slaveALDR(j,i)./nanmean(normL23slaveData{i}(j,[1:bfIndex(1)-2]));
        elseif normL23slaveData{i}(j,bfIndex(1)-1)>= normL23slaveData{i}(j,bfIndex(1)+1)
            L23slaveALDR(j,i) = (L23slaveALDR(j,i)+normL23slaveData{i}(j,bfIndex(1)-1))./2;
            L23slaveALDR(j,i) = L23slaveALDR(j,i)./nanmean(normL23slaveData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL23slaveData{i},2)]));
        elseif normL23slaveData{i}(j,bfIndex(1)-1)< normL23slaveData{i}(j,bfIndex(1)+1)
            L23slaveALDR(j,i) = (L23slaveALDR(j,i)+normL23slaveData{i}(j,bfIndex(1)+1))./2;
            L23slaveALDR(j,i) = L23slaveALDR(j,i)./nanmean(normL23slaveData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL23slaveData{i},2)]));
        end
    end
    
end
L23ALDR = ((L23masterALDR-L23slaveALDR)./L23slaveALDR)*100;
L23ALDR(L23ALDR == inf) = NaN;

for i = 1:length(L4masterData)
%     tmp = find(nanmean(L4masterData{i},1) == max(nanmean(L4masterData{i},1)));
%     BFL4masterData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L4masterData{i},1)
        normL4masterData{i}(j,:) = (L4masterData{i}(j,:)./max(L4masterData{i}(j,:)));
        bfIndex = find(normL4masterData{i}(j,:) == 1);
        L4masterALDR(j,i) = normL4masterData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L4masterALDR(j,i) = (L4masterALDR(j,i)+normL4masterData{i}(j,bfIndex(1)+1))./2;
            L4masterALDR(j,i) = L4masterALDR(j,i)./nanmean(normL4masterData{i}(j,[bfIndex(1)+2:size(normL4masterData{i},2)]));
        elseif (bfIndex(1)+1)> size(L4masterData{i},2)
            L4masterALDR(j,i) = (L4masterALDR(j,i)+normL4masterData{i}(j,bfIndex(1)-1))./2;
            L4masterALDR(j,i) = L4masterALDR(j,i)./nanmean(normL4masterData{i}(j,[1:bfIndex(1)-2]));
        elseif normL4masterData{i}(j,bfIndex(1)-1)>= normL4masterData{i}(j,bfIndex(1)+1)
            L4masterALDR(j,i) = (L4masterALDR(j,i)+normL4masterData{i}(j,bfIndex(1)-1))./2;
            L4masterALDR(j,i) = L4masterALDR(j,i)./nanmean(normL4masterData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL4masterData{i},2)]));
        elseif normL4masterData{i}(j,bfIndex(1)-1)< normL4masterData{i}(j,bfIndex(1)+1)
            L4masterALDR(j,i) = (L4masterALDR(j,i)+normL4masterData{i}(j,bfIndex(1)+1))./2;
            L4masterALDR(j,i) = L4masterALDR(j,i)./nanmean(normL4masterData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL4masterData{i},2)]));
        end
    end
    
end
for i = 1:length(L4slaveData)
%     tmp = find(nanmean(L4slaveData{i},1) == max(nanmean(L4slaveData{i},1)));
%     BFL4slaveData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L4slaveData{i},1)
        normL4slaveData{i}(j,:) = (L4slaveData{i}(j,:)./max(L4slaveData{i}(j,:)));
        bfIndex = find(normL4slaveData{i}(j,:) == 1);
        L4slaveALDR(j,i) = normL4slaveData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L4slaveALDR(j,i) = (L4slaveALDR(j,i)+normL4slaveData{i}(j,bfIndex(1)+1))./2;
            L4slaveALDR(j,i) = L4slaveALDR(j,i)./nanmean(normL4slaveData{i}(j,[bfIndex(1)+2:size(normL4slaveData{i},2)]));
        elseif (bfIndex(1)+1)> size(L4slaveData{i},2)
            L4slaveALDR(j,i) = (L4slaveALDR(j,i)+normL4slaveData{i}(j,bfIndex(1)-1))./2;
            L4slaveALDR(j,i) = L4slaveALDR(j,i)./nanmean(normL4slaveData{i}(j,[1:bfIndex(1)-2]));
        elseif normL4slaveData{i}(j,bfIndex(1)-1)>= normL4slaveData{i}(j,bfIndex(1)+1)
            L4slaveALDR(j,i) = (L4slaveALDR(j,i)+normL4slaveData{i}(j,bfIndex(1)-1))./2;
            L4slaveALDR(j,i) = L4slaveALDR(j,i)./nanmean(normL4slaveData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL4slaveData{i},2)]));
        elseif normL4slaveData{i}(j,bfIndex(1)-1)< normL4slaveData{i}(j,bfIndex(1)+1)
            L4slaveALDR(j,i) = (L4slaveALDR(j,i)+normL4slaveData{i}(j,bfIndex(1)+1))./2;
            L4slaveALDR(j,i) = L4slaveALDR(j,i)./nanmean(normL4slaveData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL4slaveData{i},2)]));
        end
    end
    
end
L4ALDR = ((L4masterALDR-L4slaveALDR)./L4slaveALDR)*100;
L4ALDR(L4ALDR == inf) = NaN;

for i = 1:length(L5masterData)
%     tmp = find(nanmean(L5masterData{i},1) == max(nanmean(L5masterData{i},1)));
%     BFL5masterData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L5masterData{i},1)
        normL5masterData{i}(j,:) = (L5masterData{i}(j,:)./max(L5masterData{i}(j,:)));
        bfIndex = find(normL5masterData{i}(j,:) == 1);
        L5masterALDR(j,i) = normL5masterData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L5masterALDR(j,i) = (L5masterALDR(j,i)+normL5masterData{i}(j,bfIndex(1)+1))./2;
            L5masterALDR(j,i) = L5masterALDR(j,i)./nanmean(normL5masterData{i}(j,[bfIndex(1)+2:size(normL5masterData{i},2)]));
        elseif (bfIndex(1)+1)> size(L5masterData{i},2)
            L5masterALDR(j,i) = (L5masterALDR(j,i)+normL5masterData{i}(j,bfIndex(1)-1))./2;
            L5masterALDR(j,i) = L5masterALDR(j,i)./nanmean(normL5masterData{i}(j,[1:bfIndex(1)-2]));
        elseif normL5masterData{i}(j,bfIndex(1)-1)>= normL5masterData{i}(j,bfIndex(1)+1)
            L5masterALDR(j,i) = (L5masterALDR(j,i)+normL5masterData{i}(j,bfIndex(1)-1))./2;
            L5masterALDR(j,i) = L5masterALDR(j,i)./nanmean(normL5masterData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL5masterData{i},2)]));
        elseif normL5masterData{i}(j,bfIndex(1)-1)< normL5masterData{i}(j,bfIndex(1)+1)
            L5masterALDR(j,i) = (L5masterALDR(j,i)+normL5masterData{i}(j,bfIndex(1)+1))./2;
            L5masterALDR(j,i) = L5masterALDR(j,i)./nanmean(normL5masterData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL5masterData{i},2)]));
        end
    end
    
end
for i = 1:length(L5slaveData)
%     tmp = find(nanmean(L5slaveData{i},1) == max(nanmean(L5slaveData{i},1)));
%     BFL5slaveData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L5slaveData{i},1)
        normL5slaveData{i}(j,:) = (L5slaveData{i}(j,:)./max(L5slaveData{i}(j,:)));
        bfIndex = find(normL5slaveData{i}(j,:) == 1);
        L5slaveALDR(j,i) = normL5slaveData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L5slaveALDR(j,i) = (L5slaveALDR(j,i)+normL5slaveData{i}(j,bfIndex(1)+1))./2;
            L5slaveALDR(j,i) = L5slaveALDR(j,i)./nanmean(normL5slaveData{i}(j,[bfIndex(1)+2:size(normL5slaveData{i},2)]));
        elseif (bfIndex(1)+1)> size(L5slaveData{i},2)
            L5slaveALDR(j,i) = (L5slaveALDR(j,i)+normL5slaveData{i}(j,bfIndex(1)-1))./2;
            L5slaveALDR(j,i) = L5slaveALDR(j,i)./nanmean(normL5slaveData{i}(j,[1:bfIndex(1)-2]));
        elseif normL5slaveData{i}(j,bfIndex(1)-1)>= normL5slaveData{i}(j,bfIndex(1)+1)
            L5slaveALDR(j,i) = (L5slaveALDR(j,i)+normL5slaveData{i}(j,bfIndex(1)-1))./2;
            L5slaveALDR(j,i) = L5slaveALDR(j,i)./nanmean(normL5slaveData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL5slaveData{i},2)]));
        elseif normL5slaveData{i}(j,bfIndex(1)-1)< normL5slaveData{i}(j,bfIndex(1)+1)
            L5slaveALDR(j,i) = (L5slaveALDR(j,i)+normL5slaveData{i}(j,bfIndex(1)+1))./2;
            L5slaveALDR(j,i) = L5slaveALDR(j,i)./nanmean(normL5slaveData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL5slaveData{i},2)]));
        end
    end
    
end
L5ALDR = ((L5masterALDR-L5slaveALDR)./L5slaveALDR)*100;
L5ALDR(L5ALDR == inf) = NaN;

for i = 1:length(L6masterData)
%     tmp = find(nanmean(L6masterData{i},1) == max(nanmean(L6masterData{i},1)));
%     BFL6masterData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L6masterData{i},1)
        normL6masterData{i}(j,:) = (L6masterData{i}(j,:)./max(L6masterData{i}(j,:)));
        bfIndex = find(normL6masterData{i}(j,:) == 1);
        L6masterALDR(j,i) = normL6masterData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L6masterALDR(j,i) = (L6masterALDR(j,i)+normL6masterData{i}(j,bfIndex(1)+1))./2;
            L6masterALDR(j,i) = L6masterALDR(j,i)./nanmean(normL6masterData{i}(j,[bfIndex(1)+2:size(normL6masterData{i},2)]));
        elseif (bfIndex(1)+1)> size(L6masterData{i},2)
            L6masterALDR(j,i) = (L6masterALDR(j,i)+normL6masterData{i}(j,bfIndex(1)-1))./2;
            L6masterALDR(j,i) = L6masterALDR(j,i)./nanmean(normL6masterData{i}(j,[1:bfIndex(1)-2]));
        elseif normL6masterData{i}(j,bfIndex(1)-1)>= normL6masterData{i}(j,bfIndex(1)+1)
            L6masterALDR(j,i) = (L6masterALDR(j,i)+normL6masterData{i}(j,bfIndex(1)-1))./2;
            L6masterALDR(j,i) = L6masterALDR(j,i)./nanmean(normL6masterData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL6masterData{i},2)]));
        elseif normL6masterData{i}(j,bfIndex(1)-1)< normL6masterData{i}(j,bfIndex(1)+1)
            L6masterALDR(j,i) = (L6masterALDR(j,i)+normL6masterData{i}(j,bfIndex(1)+1))./2;
            L6masterALDR(j,i) = L6masterALDR(j,i)./nanmean(normL6masterData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL6masterData{i},2)]));
        end
    end
    
end
for i = 1:length(L6slaveData)
%     tmp = find(nanmean(L6slaveData{i},1) == max(nanmean(L6slaveData{i},1)));
%     BFL6slaveData(i) = tmp(1);
    clear tmp;
    for j = 1:size(L6slaveData{i},1)
        normL6slaveData{i}(j,:) = (L6slaveData{i}(j,:)./max(L6slaveData{i}(j,:)));
        bfIndex = find(normL6slaveData{i}(j,:) == 1);
        L6slaveALDR(j,i) = normL6slaveData{i}(j,bfIndex(1));
        if (bfIndex(1)-1)== 0
            L6slaveALDR(j,i) = (L6slaveALDR(j,i)+normL6slaveData{i}(j,bfIndex(1)+1))./2;
            L6slaveALDR(j,i) = L6slaveALDR(j,i)./nanmean(normL6slaveData{i}(j,[bfIndex(1)+2:size(normL6slaveData{i},2)]));
        elseif (bfIndex(1)+1)> size(L6slaveData{i},2)
            L6slaveALDR(j,i) = (L6slaveALDR(j,i)+normL6slaveData{i}(j,bfIndex(1)-1))./2;
            L6slaveALDR(j,i) = L6slaveALDR(j,i)./nanmean(normL6slaveData{i}(j,[1:bfIndex(1)-2]));
        elseif normL6slaveData{i}(j,bfIndex(1)-1)>= normL6slaveData{i}(j,bfIndex(1)+1)
            L6slaveALDR(j,i) = (L6slaveALDR(j,i)+normL6slaveData{i}(j,bfIndex(1)-1))./2;
            L6slaveALDR(j,i) = L6slaveALDR(j,i)./nanmean(normL6slaveData{i}(j,[1:bfIndex(1)-2 bfIndex(1)+1:size(normL6slaveData{i},2)]));
        elseif normL6slaveData{i}(j,bfIndex(1)-1)< normL6slaveData{i}(j,bfIndex(1)+1)
            L6slaveALDR(j,i) = (L6slaveALDR(j,i)+normL6slaveData{i}(j,bfIndex(1)+1))./2;
            L6slaveALDR(j,i) = L6slaveALDR(j,i)./nanmean(normL6slaveData{i}(j,[1:bfIndex(1)-1 bfIndex(1)+2:size(normL6slaveData{i},2)]));
        end
    end
    
end
L6ALDR = ((L6masterALDR-L6slaveALDR)./L6slaveALDR)*100;
L6ALDR(L6ALDR == inf) = NaN;





%% Mean and SEM

    meanL23gainData = fliplr(mean(L23gainData,2)');
    semL23gainData = fliplr(std(L23gainData')./sqrt(size(L23gainData,2))');
    meanL23ALDR = fliplr(nanmean(L23ALDR,2)');
    semL23ALDR = fliplr((nanstd(L23ALDR'))./sqrt(size(L23ALDR,2)'));
    meanL23psth = mean(L23psth);
    semL23psth = std(L23psth)./sqrt(size(L23psth,1));
    
    meanL4gainData = fliplr(mean(L4gainData,2)');
    semL4gainData = fliplr(std(L4gainData')./sqrt(size(L4gainData,2))');
    meanL4ALDR = fliplr(nanmean(L4ALDR,2)');
    semL4ALDR = fliplr((nanstd(L4ALDR'))./sqrt(size(L4ALDR,2)'));
    meanL4psth = mean(L4psth);
    semL4psth = std(L4psth)./sqrt(size(L4psth,1));
     
    meanL5gainData = fliplr(mean(L5gainData,2)');
    semL5gainData = fliplr(std(L5gainData')./sqrt(size(L5gainData,2))');
    meanL5ALDR = fliplr(nanmean(L5ALDR,2)');
    semL5ALDR = fliplr((nanstd(L5ALDR'))./sqrt(size(L5ALDR,2)'));
    meanL5psth = mean(L5psth);
    semL5psth = std(L5psth)./sqrt(size(L5psth,1));
    
    meanL6gainData = fliplr(mean(L6gainData,2)');
    semL6gainData = fliplr(std(L6gainData')./sqrt(size(L6gainData,2))');
    meanL6ALDR = fliplr(nanmean(L6ALDR,2)');
    semL6ALDR = fliplr((nanstd(L6ALDR'))./sqrt(size(L6ALDR,2)'));
    meanL6psth = mean(L6psth);
    semL6psth = std(L6psth)./sqrt(size(L6psth,1));
    

% Together
figure
alpha = .9;

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
area([0:25:400,400:-25:0],[meanL23ALDR-semL23ALDR,meanL23ALDR(end:-1:1)+semL23ALDR(end:-1:1)],'facecolor',[ 0 alpha 0],'linestyle','none');
hold on
area([0:25:400,400:-25:0],[meanL23gainData-semL23gainData,meanL23gainData(end:-1:1)+semL23gainData(end:-1:1)],'facecolor',[ alpha 0 0],'linestyle','none');
area([0:1:400,400:-1:0],[meanL23psth(1:401)-semL23psth(1:401),meanL23psth(401:-1:1)+semL23psth(401:-1:1)]*550,'facecolor',[0 0 1],'linestyle','none');
% plot(0:50:800,fliplr(L23gainData')','color',[0 0 0]+alpha)
plot(0:25:400,meanL23ALDR,'linewidth',3,'color',[0 .3 0])
plot(0:25:400,meanL23gainData,'linewidth',3,'color',[.3 0 0])
plot(0:400,meanL23psth(1:401)*550,'linewidth',3,'color',[0 0 .3])
plot(0:25:400,zeros(17,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([0 400 -120 650])
%plot(0:25:400,fliplr(gSelectivity(:,1)')*100)

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,2)
area([0:25:400,400:-25:0],[meanL4ALDR-semL4ALDR,meanL4ALDR(end:-1:1)+semL4ALDR(end:-1:1)],'facecolor',[ 0 alpha 0],'linestyle','none');
hold on
area([0:25:400,400:-25:0],[meanL4gainData-semL4gainData,meanL4gainData(end:-1:1)+semL4gainData(end:-1:1)],'facecolor',[ alpha 0 0],'linestyle','none');
area([0:1:400,400:-1:0],[meanL4psth(1:401)-semL4psth(1:401),meanL4psth(401:-1:1)+semL4psth(401:-1:1)]*550,'facecolor',[0 0 1],'linestyle','none');
% plot(0:50:800,fliplr(L4gainData')','color',[0 0 0]+alpha)
plot(0:25:400,meanL4ALDR,'linewidth',3,'color',[0 .3 0])
plot(0:25:400,meanL4gainData,'linewidth',3,'color',[.3 0 0])
plot(0:400,meanL4psth(1:401)*550,'linewidth',3,'color',[0 0 .3])
plot(0:25:400,zeros(17,1),'k--')
title('Layer 4')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([0 400 -120 650])
%plot(0:25:400,fliplr(gSelectivity(:,2)')*100)

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,3)
area([0:25:400,400:-25:0],[meanL5ALDR-semL5ALDR,meanL5ALDR(end:-1:1)+semL5ALDR(end:-1:1)],'facecolor',[ 0 alpha 0],'linestyle','none');
hold on
area([0:25:400,400:-25:0],[meanL5gainData-semL5gainData,meanL5gainData(end:-1:1)+semL5gainData(end:-1:1)],'facecolor',[ alpha 0 0],'linestyle','none');
area([0:1:400,400:-1:0],[meanL5psth(1:401)-semL5psth(1:401),meanL5psth(401:-1:1)+semL5psth(401:-1:1)]*550,'facecolor',[0 0 1],'linestyle','none');
% plot(0:50:800,fliplr(L5gainData')','color',[0 0 0]+alpha)
plot(0:25:400,meanL5ALDR,'linewidth',3,'color',[0 .3 0])
plot(0:25:400,meanL5gainData,'linewidth',3,'color',[.3 0 0])
plot(0:400,meanL5psth(1:401)*550,'linewidth',3,'color',[0 0 .3])
plot(0:25:400,zeros(17,1),'k--')
title('Layer 5')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([0 400 -120 650])
%plot(0:25:400,fliplr(gSelectivity(:,3)')*100)

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,4)
area([0:25:400,400:-25:0],[meanL6ALDR-semL6ALDR,meanL6ALDR(end:-1:1)+semL6ALDR(end:-1:1)],'facecolor',[ 0 alpha 0],'linestyle','none');
hold on
area([0:25:400,400:-25:0],[meanL6gainData-semL6gainData,meanL6gainData(end:-1:1)+semL6gainData(end:-1:1)],'facecolor',[ alpha 0 0],'linestyle','none');
area([0:1:400,400:-1:0],[meanL6psth(1:401)-semL6psth(1:401),meanL6psth(401:-1:1)+semL6psth(401:-1:1)]*550,'facecolor',[0 0 1],'linestyle','none');
% plot(0:50:800,fliplr(L6gainData')','color',[0 0 0]+alpha)
plot(0:25:400,meanL6ALDR,'linewidth',3,'color',[0 .3 0])
plot(0:25:400,meanL6gainData,'linewidth',3,'color',[.3 0 0])
plot(0:400,meanL6psth(1:401)*550,'linewidth',3,'color',[0 0 .3])
plot(0:25:400,zeros(17,1),'k--')
title('Layer 6')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([0 400 -120 650])   
%plot(0:25:400,fliplr(gSelectivity(:,4)')*100)
%% Plots
%{

figure
rgbVals = 0:0.0625:1;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.03], [0.05 0.06], [0.05 0.02]);

subplot(2,2,1)
h = 1;
for i = 1:size(L23gainData,1)
    plot(meanToneL23slaveData(i,:), meanToneL23masterData(i,:),'*','color', [rgbVals(h),0,1-rgbVals(h)])
    h = h+1;
    hold on
end
hold on
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 2/3')
xlabel('Tone alone')
ylabel('Tone + laser')


subplot(2,2,2)
h = 1;
for i = 1:size(L4gainData,1)
    plot(meanToneL4slaveData(i,:), meanToneL4masterData(i,:),'*','color', [rgbVals(h),0,1-rgbVals(h)])
    h = h+1;
    hold on
end
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 4')
xlabel('Tone alone')
ylabel('Tone + laser')

subplot(2,2,5)
h = 1;
for i = 1:size(L5gainData,1)
    plot(meanToneL5slaveData(i,:), meanToneL5masterData(i,:),'*','color', [rgbVals(h),0,1-rgbVals(h)])
    h = h+1;
    hold on
end
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 5')
xlabel('Tone alone')
ylabel('Tone + laser')

subplot(2,2,6)
h = 1;
for i = 1:size(L6gainData,1)
    plot(meanToneL6slaveData(i,:), meanToneL6masterData(i,:),'*','color', [rgbVals(h),0,1-rgbVals(h)])
    h = h+1;
    hold on
end
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 6')
xlabel('Tone alone')
ylabel('Tone + laser')


figure
alpha = .8;

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
area([0:50:800,800:-50:0],[meanL23gainData-semL23gainData,meanL23gainData(end:-1:1)+semL23gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L23gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL23gainData,'linewidth',3,'color',[0 0 0])
%plot(-100:800,meanL23psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain (% increase)')
set(gca,'box','off')
axis([0 800 -100 250])

subplot(2,2,2)
area([0:50:800,800:-50:0],[meanL4gainData-semL4gainData,meanL4gainData(end:-1:1)+semL4gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L4gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL4gainData,'linewidth',3,'color',[0 0 0])
%plot(-100:800,meanL4psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 4')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain (% increase)')
set(gca,'box','off')
axis([0 800 -100 250])

subplot(2,2,3)
area([0:50:800,800:-50:0],[meanL5gainData-semL5gainData,meanL5gainData(end:-1:1)+semL5gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L5gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL5gainData,'linewidth',3,'color',[0 0 0])
%plot(-100:800,meanL5psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 5')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain (% increase)')
set(gca,'box','off')
axis([0 800 -100 250])

subplot(2,2,4)
area([0:50:800,800:-50:0],[meanL6gainData-semL6gainData,meanL6gainData(end:-1:1)+semL6gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L6gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL6gainData,'linewidth',3,'color',[0 0 0])
%plot(-100:800,meanL6psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 6')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain (% increase)')
set(gca,'box','off')
axis([0 800 -100 250])
   

% gain
figure
alpha = .8;

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
area([0:50:800,800:-50:0],[meanL23gainData-semL23gainData,meanL23gainData(end:-1:1)+semL23gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L23gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL23gainData,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL23psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 250])

subplot(2,2,2)
area([0:50:800,800:-50:0],[meanL4gainData-semL4gainData,meanL4gainData(end:-1:1)+semL4gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L4gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL4gainData,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL4psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 4')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 250])

subplot(2,2,3)
area([0:50:800,800:-50:0],[meanL5gainData-semL5gainData,meanL5gainData(end:-1:1)+semL5gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L5gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL5gainData,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL5psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 5')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain)')
set(gca,'box','off')
axis([-100 800 -100 250])

subplot(2,2,4)
area([0:50:800,800:-50:0],[meanL6gainData-semL6gainData,meanL6gainData(end:-1:1)+semL6gainData(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L6gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL6gainData,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL6psth(100:1000)*250,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 6')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 250])

% Selectivity
figure
alpha = .8;

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
area([0:50:800,800:-50:0],[meanL23ALDR-semL23ALDR,meanL23ALDR(end:-1:1)+semL23ALDR(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L23gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL23ALDR,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL23psth(100:1000)*300,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 300])

subplot(2,2,2)
area([0:50:800,800:-50:0],[meanL4ALDR-semL4ALDR,meanL4ALDR(end:-1:1)+semL4ALDR(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L4gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL4ALDR,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL4psth(100:1000)*300,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 300])

subplot(2,2,3)
area([0:50:800,800:-50:0],[meanL5ALDR-semL5ALDR,meanL5ALDR(end:-1:1)+semL5ALDR(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L5gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL5ALDR,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL5psth(100:1000)*300,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 300])

subplot(2,2,4)
area([0:50:800,800:-50:0],[meanL6ALDR-semL6ALDR,meanL6ALDR(end:-1:1)+semL6ALDR(end:-1:1)],'facecolor',[ 0 0 0]+alpha,'linestyle','none');
% plot(0:50:800,fliplr(L6gainData')','color',[0 0 0]+alpha)
hold on
plot(0:50:800,meanL6ALDR,'linewidth',3,'color',[0 0 0])
plot(-100:800,meanL6psth(100:1000)*300,'linewidth',2)
plot(-100:50:800,zeros(19,1),'k--')
title('Layer 2/3')
xlabel('Tone onset (relative to laser onset)')
ylabel('Gain')
set(gca,'box','off')
axis([-100 800 -100 300])

%}

% figure
% 
% subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);
% subplot(2,2,1)
% plot( meanToneL23slaveData, meanToneL23masterData,'*')
% hold on
% y = get(gca,'YLim')
% x = get(gca,'XLim')
% maxAx = max([x(2) y(2)])
% plot(1:maxAx,1:maxAx,'k--')
% title('Layer 2/3')
% xlabel('Tone alone')
% ylabel('Tone + laser')
% 
% subplot(2,2,2)
% plot( meanToneL4slaveData, meanToneL4masterData,'*')
% hold on
% y = get(gca,'YLim')
% x = get(gca,'XLim')
% maxAx = max([x(2) y(2)])
% plot(1:maxAx,1:maxAx,'k--')
% title('Layer 4')
% xlabel('Tone alone')
% ylabel('Tone + laser')
% 
% subplot(2,2,3)
% plot( meanToneL5slaveData, meanToneL5masterData,'*')
% hold on
% y = get(gca,'YLim')
% x = get(gca,'XLim')
% maxAx = max([x(2) y(2)])
% plot(1:maxAx,1:maxAx,'k--')
% title('Layer 5')
% xlabel('Tone alone')
% ylabel('Tone + laser')
% 
% subplot(2,2,4)
% plot( meanToneL6slaveData, meanToneL6masterData,'*')
% hold on
% y = get(gca,'YLim')
% x = get(gca,'XLim')
% maxAx = max([x(2) y(2)])
% plot(1:maxAx,1:maxAx,'k--')
% title('Layer 6')
% xlabel('Tone alone')
% ylabel('Tone + laser')
% 
% 
%     figure
% 
%     subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);
%     subplot(2,2,1)
%     for i =1:length(L23masterData)
%         plot(mean(L23slaveData{i},2),mean(L23masterData{i},2),'k*')
%         hold on
%     end
%     y = get(gca,'YLim')
%     x = get(gca,'XLim')
%     maxAx = max([x(2) y(2)])
%     plot(1:maxAx,1:maxAx,'k--')
%     title('Layer 2/3')
%     xlabel('Tone alone')
%     ylabel('Tone + laser')
% 
%     subplot(2,2,2)
%     for i =1:length(L4masterData)
%         plot(mean(L4slaveData{i},2),mean(L4masterData{i},2),'k*')
%         hold on
%     end
%     y = get(gca,'YLim')
%     x = get(gca,'XLim')
%     maxAx = max([x(2) y(2)])
%     plot(1:maxAx,1:maxAx,'k--')
%     title('Layer 4')
%     xlabel('Tone alone')
%     ylabel('Tone + laser')
% 
%     subplot(2,2,3)
%     for i =1:length(L5masterData)
%         plot(mean(L5slaveData{i},2),mean(L5masterData{i},2),'k*')
%         hold on
%     end
%     y = get(gca,'YLim')
%     x = get(gca,'XLim')
%     maxAx = max([x(2) y(2)])
%     plot(1:maxAx,1:maxAx,'k--')
%     title('Layer 5')
%     xlabel('Tone alone')
%     ylabel('Tone + laser')
% 
%     subplot(2,2,4)
%     for i =1:length(L6masterData)
%         plot(mean(L6slaveData{i},2),mean(L6masterData{i},2),'k*')
%         hold on
%     end
%     y = get(gca,'YLim')
%     x = get(gca,'XLim')
%     maxAx = max([x(2) y(2)])
%     plot(1:maxAx,1:maxAx,'k--')
%     title('Layer 6')
%     xlabel('Tone alone')
%     ylabel('Tone + laser')