for i = 1:length(SCL)
    if isempty(SCL(i).t)
        break
    end
    lastTrial = i;
end

wholeProbeTimes = [];

for i =1:lastTrial
    wholeProbeTimes = [wholeProbeTimes;SCL(i).t];
end

wholeProbeRaster = zeros((dacInt.SampleRate/10) * (wholeProbeTimes(lastTrial)/1000),1);
wholeProbeRaster(ceil(wholeProbeTimes*(dacInt.SampleRate/10000))) = 1;
wholeProbePSTH = downsample(quickPSTH(wholeProbeRaster,200),10);

figure
plot(wholeProbePSTH,'k')

9.962e5

99,6200

274

  
 for i = [7,13,4,12,5,15,2,16,1,14,3,11,6]
     clims = [min(min([masterData{i}-masterStdData{i},slaveData{i}-slaveStdData{i}])),...
         max(max([masterData{i}+masterStdData{i},slaveData{i}+slaveStdData{i}]))];
     figure
     subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.04], [0.05 0.02], [0.05 0.02]);
     for j = 1:length(masterData{i})
         subplot(9,2,j)
         errorbar(masterData{i}(j,:),masterStdData{i}(j,:),'linewidth',2)
         hold on
         errorbar(slaveData{i}(j,:),slaveStdData{i}(j,:),'color',[0,.4,.2],'linewidth',2)
         axis([1 8 clims(1) clims(2)])
         set(gca,'XTickLabel',{'4','5.7','8','11.3','16','22.6','32','45.3'});
         set(gca,'box','off')
         ylabel([num2str((j*50)+150) 'ms delay'])
     end
     % text(7.8,165,['Channel #' num2str(i)])
     subplot(9,2,1)
     title(['Channel #' num2str(i)])
     subplot(9,2,16)
     xlabel('Frequency (KHz)')
     subplot(9,2,17)
     xlabel('Frequency (KHz)')
     
 end