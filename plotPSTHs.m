channel = 3;
% Pick a channel and average the PSTHs across repetitions
repPSTHs = repPsthData{channel};
for i = 1:17

    repAvgPSTHs{i} = (repPSTHs{i,1}+repPSTHs{i,2}+repPSTHs{i,3}+repPSTHs{i,4}+repPSTHs{i,5}+repPSTHs{i,6}+repPSTHs{i,7}+repPSTHs{i,8})/8;

end

% Plot PSTHs averaged above and in the firstAnalysis code to make sure they
% converge
 figure
 spacer = 1.5;
for i =1:17
    plot((squeeze(mean(psthData{channel}(i,:,:),2))*1000)+(i*spacer),'k');
    plot(mean(repAvgPSTHs{i}*1000)+(i*spacer),'color',[0 0 1])
    hold on
    h = patch([(150+(i*50)) (150+(i*50)) (550+(i*50)) (550+(i*50))],[(.4+(i*spacer)) (2.4+(i*spacer)) (2.4+(i*spacer)) (.4+(i*spacer))],'b');
    set(h, 'EdgeColor','none')
    set(h, 'FaceAlpha',.3)
    hold on
end

% Plot mean PSTHs with individual repetitions in gray, means in black, and
% laser as light blue patch.
figure
spacer = 13;
hold on
for i =1:17
    alpha = .8;
    plot(squeeze(repAvgPSTHs{i}*1000)'+(i*spacer),'color',[0 0 0]+alpha)
    plot(mean(repAvgPSTHs{i}*1000)+(i*spacer),'linewidth',2,'color',[0 0 0])
    h = patch([(150+(i*50)) (150+(i*50)) (550+(i*50)) (550+(i*50))],[(.4+(i*spacer)) (2.4+(i*spacer)) (2.4+(i*spacer)) (.4+(i*spacer))],'b');
    set(h, 'EdgeColor','none')
    set(h, 'FaceAlpha',.3)
    plot(1000*ones(251,1),0:250,'--','linewidth',2)
end
title(['Channel ' num2str(channel)])
hold off

