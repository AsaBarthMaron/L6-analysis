%function[L23masterData, L4masterData, L5masterData, L6masterData] = quantifyGain(indices,toneResponses,startingBins,channelAssignment)

[masterDataExp,slaveDataExp] = laserDelay_tuningCurve_secondAnalysis(indices,toneResponses,startingBins,0)
close all
L23masterData = [];
L4masterData = [];
L5masterData = [];
L6masterData = [];
L23slaveData = [];
L4slaveData = [];
L5slaveData = [];
L6slaveData = [];
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
        elseif (max(channelAssignment(k,4:5) == ch)) && (max(toneResponses{k} == ch))
            L4masterData = [L4masterData masterDataExp(k,ch)];
            L4slaveData = [L4slaveData slaveDataExp(k,ch)];
        elseif (max(channelAssignment(k,6:8) == ch)) && (max(toneResponses{k} == ch))
            L5masterData = [L5masterData masterDataExp(k,ch)];
            L5slaveData = [L5slaveData slaveDataExp(k,ch)];
        elseif (max(channelAssignment(k,9:11) == ch)) && (max(toneResponses{k} == ch))
            L6masterData = [L6masterData masterDataExp(k,ch)];
            L6slaveData = [L6slaveData slaveDataExp(k,ch)];
        end
    end
end

for i =1:length(L23masterData)
    L23gainData(:,i) = (mean(L23slaveData{i},2)-mean(L23masterData{i},2))...
                     ./(mean(L23slaveData{i},2)+mean(L23masterData{i},2));
    meanToneL23masterData(:,i) = mean(L23masterData{i},2);
    meanToneL23slaveData(:,i) = mean(L23slaveData{i},2);
end

for i =1:length(L4masterData)
    L4gainData(:,i) = (mean(L4slaveData{i},2)-mean(L4masterData{i},2))...
                    ./(mean(L4slaveData{i},2)+mean(L4masterData{i},2));
    meanToneL4masterData(:,i) = mean(L4masterData{i},2);
    meanToneL4slaveData(:,i) = mean(L4slaveData{i},2);
end

for i =1:length(L5masterData)
    L5gainData(:,i) = (mean(L5slaveData{i},2)-mean(L5masterData{i},2))...
                    ./(mean(L5slaveData{i},2)+mean(L5masterData{i},2));
    meanToneL5masterData(:,i) = mean(L5masterData{i},2);
    meanToneL5slaveData(:,i) = mean(L5slaveData{i},2);
end

for i =1:length(L6masterData)
    L6gainData(:,i) = (mean(L6slaveData{i},2)-mean(L6masterData{i},2))...
                    ./(mean(L6slaveData{i},2)+mean(L6masterData{i},2));
    meanToneL6masterData(:,i) = mean(L6masterData{i},2);
    meanToneL6slaveData(:,i) = mean(L6slaveData{i},2);
end


figure
rgbVals = 0:0.0625:1;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);

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

subplot(2,2,3)
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

subplot(2,2,4)
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

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
plot(-800:50:0,L23gainData)
hold on
plot(-800:50:0,zeros(17,1),'k--')
title('Layer 2/3')
xlabel('Laser onset (relative to tone)')
ylabel('Gain')

subplot(2,2,2)
plot(-800:50:0,L4gainData)
hold on
plot(-800:50:0,zeros(17,1),'k--')
title('Layer 4')
xlabel('Laser onset (relative to tone)')
ylabel('Gain')

subplot(2,2,3)
plot(-800:50:0,L5gainData)
hold on
plot(-800:50:0,zeros(17,1),'k--')
title('Layer 5')
xlabel('Laser onset (relative to tone)')
ylabel('Gain')

subplot(2,2,4)
plot(-800:50:0,L6gainData)
hold on
plot(-800:50:0,zeros(17,1),'k--')
title('Layer 6')
xlabel('Laser onset (relative to tone)')
ylabel('Gain')


figure

subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);
subplot(2,2,1)
plot( meanToneL23slaveData, meanToneL23masterData,'*')
hold on
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 2/3')
xlabel('Tone alone')
ylabel('Tone + laser')

subplot(2,2,2)
plot( meanToneL4slaveData, meanToneL4masterData,'*')
hold on
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 4')
xlabel('Tone alone')
ylabel('Tone + laser')

subplot(2,2,3)
plot( meanToneL5slaveData, meanToneL5masterData,'*')
hold on
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 5')
xlabel('Tone alone')
ylabel('Tone + laser')

subplot(2,2,4)
plot( meanToneL6slaveData, meanToneL6masterData,'*')
hold on
y = get(gca,'YLim')
x = get(gca,'XLim')
maxAx = max([x(2) y(2)])
plot(1:maxAx,1:maxAx,'k--')
title('Layer 6')
xlabel('Tone alone')
ylabel('Tone + laser')


    figure

    subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.06], [0.05 0.02]);
    subplot(2,2,1)
    for i =1:length(L23masterData)
        plot(mean(L23slaveData{i},2),mean(L23masterData{i},2),'k*')
        hold on
    end
    y = get(gca,'YLim')
    x = get(gca,'XLim')
    maxAx = max([x(2) y(2)])
    plot(1:maxAx,1:maxAx,'k--')
    title('Layer 2/3')
    xlabel('Tone alone')
    ylabel('Tone + laser')

    subplot(2,2,2)
    for i =1:length(L4masterData)
        plot(mean(L4slaveData{i},2),mean(L4masterData{i},2),'k*')
        hold on
    end
    y = get(gca,'YLim')
    x = get(gca,'XLim')
    maxAx = max([x(2) y(2)])
    plot(1:maxAx,1:maxAx,'k--')
    title('Layer 4')
    xlabel('Tone alone')
    ylabel('Tone + laser')

    subplot(2,2,3)
    for i =1:length(L5masterData)
        plot(mean(L5slaveData{i},2),mean(L5masterData{i},2),'k*')
        hold on
    end
    y = get(gca,'YLim')
    x = get(gca,'XLim')
    maxAx = max([x(2) y(2)])
    plot(1:maxAx,1:maxAx,'k--')
    title('Layer 5')
    xlabel('Tone alone')
    ylabel('Tone + laser')

    subplot(2,2,4)
    for i =1:length(L6masterData)
        plot(mean(L6slaveData{i},2),mean(L6masterData{i},2),'k*')
        hold on
    end
    y = get(gca,'YLim')
    x = get(gca,'XLim')
    maxAx = max([x(2) y(2)])
    plot(1:maxAx,1:maxAx,'k--')
    title('Layer 6')
    xlabel('Tone alone')
    ylabel('Tone + laser')